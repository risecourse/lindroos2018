#
'''
MSN model used in Lindroos et al., (2018). Frontiers

Robert Lindroos (RL) <robert.lindroos at ki.se>
 
The MSN class and most channels were implemented by 
Alexander Kozlov <akozlov at kth.se>
with updates by RL

Implemented in colaboration with Kai Du <kai.du at ki.se>
'''


from __future__ import print_function, division
from neuron import h
#import mpi4py as MPI
import numpy                as np
import plot_functions       as fun
import MSN_builder          as build
import json


h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')



# global result dict
RES = {}


# initial and maximal values of substrates seen with a DA transient of 500 nA amp
# and 500 ms time constant.
with open('./substrates.json') as file:
    SUBSTRATES = json.load(file)


    
def alpha(tstart, gmax, tau):
    # calc and returns a "magnitude" using an alpha function -> used for [DA] in cascade
    
    v = 1 - (h.t - tstart) / tau
    e = np.exp(v)
    mag = gmax * (h.t - tstart) / tau * e
    
    return mag
    
    
    
def save_vector(t, v, outfile):
    
    with open(outfile, "w") as out:
        for time, y in zip(t, v):
            out.write("%g %g\n" % (time, y))
            
            
            
def calc_rand_Modulation(mod_list, range_list=False):
    
    '''
    uses numpy to draws random modulation factors in range [0,2],
    from a uniform distribution, for each channel in mod_list.
    
    The factors can also be linearly mapped to an arbitrary interval. 
    This is done if a range_list is given.
    
    mod_list     = list of channels to be modulated
    range_list   = list of [min, max] values to be used in modulation. 
                    Must have same length as mod_list.
    '''
    
    mod_factors = []
    
    A=0
    B=2
    
    for i,channel in enumerate(mod_list):
        
        factor = 2.0 * np.random.uniform()
    
        if range_list:
            
            a       = range_list[i][0]
            b       = range_list[i][1]
            
            factor = (b-a) / (B-A) * (factor-A) + a
       
        mod_factors.append(factor)
        
    return mod_factors 
        
         


def make_random_synapse(ns, nc, Syn, sec, x,               \
                Type='glut',                    \
                NS_start=1,                     \
                NS_interval=1000.0/18.0,        \
                NS_noise=1.0,                   \
                NS_number=1000,                 \
                S_AN_ratio=1.0,                 \
                S_tau_dep=100,                  \
                S_U=1,                          \
                S_e=-60,                        \
                S_tau1=0.25,                    \
                S_tau2=3.75,                    \
                NC_delay=1,                     \
                NC_conductance=0.6e-3,          \
                NC_threshold=0.1                ):
    
    '''
    creates a synapse in the segment closest to x in section sec.
    
    NS-arguments are used to define a NetStim object
    
    S-arguemts are used to specify the synapse mechanism
    
    NC-arguemts are used to define the NetCon object
    '''
    
    # create/set synapse in segment x of section
    if Type == 'glut':
        key                 = sec
        Syn[key]            = h.tmGlut(x, sec=sec)
        Syn[key].nmda_ratio = S_AN_ratio
        Syn[key].tauR       = S_tau_dep
        Syn[key].U          = S_U
        
    elif Type in ['expSyn2', 'tmgabaa', 'gaba']:
        
        key                 = sec.name() + '_gaba'
        
        if Type == 'expSyn2':
            Syn[key]            = h.Exp2Syn(x, sec=sec)
            Syn[key].tau1       = S_tau1
            Syn[key].tau2       = S_tau2 
        elif Type == 'tmgabaa':
            Syn[key]            = h.tmGabaA(x, sec=sec)
            
        Syn[key].e          = S_e
        
         
    
    # create NetStim object
    ns[key]             = h.NetStim()
    ns[key].start       = NS_start
    ns[key].interval    = NS_interval # mean interval between two spikes in ms
    ns[key].noise       = NS_noise
    ns[key].number      = NS_number

    
    
    # create NetCon object
    nc[key]             = h.NetCon(ns[sec],Syn[sec])
    nc[key].delay       = NC_delay
    nc[key].weight[0]   = NC_conductance
    nc[key].threshold   = NC_threshold
    



def set_rand_synapse(channel_list, base_mod, max_mod, range_list=[[0.75,1.5],[0.75,1.5]]):   
    
    syn_fact = calc_rand_Modulation(channel_list, range_list=range_list)
        
    # normalize factors to max-value of pointer substrate
    normalized_factors     = []
    for i,mech in enumerate(channel_list):
        
        normalized_factors.append( (syn_fact[i] - 1) / (max_mod - base_mod)  ) 
        
    return syn_fact, normalized_factors     
 



   
#=========================================================================================


# in the dynamimcal modulation, the channels are connected to one substrate of the cascade.
# base modulation (control) is assumed for base value of the substrate and full modulation
# is assumed when the substrate level reaches its maximal value. Linear scaling is used 
# between these points.
def main(par="./params-msn.json", \
                            sim='vm',       \
                            amp=0.265,      \
                            run=None,       \
                            dynMod=1,       \
                            simDur=7000,    \
                            stimDur=900,    \
                            target=None,    \
                            not2mod=[] ): 
    
    print('iter:', run)
    
    # initiate cell
    cell = build.MSN(params=par)
    
        
    # set cascade
    casc = h.D1_reduced_cascade2_0(0.5, sec=cell.soma)
    
    
    # set pointer (the pka proxy, Target1p, is used in this example).
    cmd         = 'pointer = casc._ref_'+target
    exec(cmd)
    base_mod    = SUBSTRATES[target][0]
    max_mod     = SUBSTRATES[target][1]
        
    
    # set edge of soma as reference for distance 
    h.distance(1, sec=h.soma[0])
    
    
    # set current injection (zero amplitude in this example)
    stim        = h.IClamp(0.5, sec=cell.soma)
    stim.amp    = amp  
    stim.delay  = 100
    stim.dur    = stimDur                
    
    
    # record vectors
    tm = h.Vector()
    tm.record(h._ref_t)
    vm = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    
    # record substrate concentrations
    if run == 0:
        pka = h.Vector()
        pka.record(casc._ref_Target1p)
        camp = h.Vector()
        camp.record(casc._ref_cAMP)
        gprot = h.Vector()
        gprot.record(casc._ref_D1RDAGolf)
        gbg   = h.Vector()
        gbg.record(casc._ref_Gbgolf)
    
    
    # parameters for DA transient
    da_peak   = 500     # concentration     [nM]
    da_tstart = 1000    # stimulation time  [ms]
    da_tau    = 500     # time constant     [ms]
    
    tstop = simDur      #                   [ms]
    # dt = default value; 0.025 ms (25 us)
    
    
    # channels to modulate channels
    mod_list    = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can' ] 
      
    
    
    
    # draw mod factors from [min, max] ranges (as percent of control). 
    # Channel ranges are in the following order:
    # ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can' ]
    mod_fact = calc_rand_Modulation(mod_list, range_list=[[0.60,0.80],  \
                                                          [0.65,0.85],  \
                                                          [0.75,0.85],  \
                                                          [0.85,1.25],  \
                                                          [1.0,2.0],    \
                                                          [1.0,2.0],    \
                                                          [0.0,1.0]]    )
    
    
    
    
    # normalize factors to target values seen in simulation by formula:
    #
    #   modulation  = 1   +      f * (substrate     - initial_substrate)
    #   f           = (factor - 1) / (max_substrate - initial_substrate)
    #
    # so that the base value correspond to no modulation and and the maximal substrate (target)
    # value corresponds to the maximal value (given by the factor).
    factors     = []
    for i,mech in enumerate(mod_list):
        
        factor  = (mod_fact[i] - 1) / (max_mod - base_mod) #2317.1
        factors.append(factor)
        
            
    
    # set pointers 
    for sec in h.allsec():
        
        for seg in sec:
            
            # naf and kas are distributed to all sections
            h.setpointer(pointer, 'pka', seg.kas )
            h.setpointer(pointer, 'pka', seg.naf )
            
            
            if sec.name().find('axon') < 0:    
                
                
                # these channels are not in the axon section
                h.setpointer(pointer, 'pka', seg.kaf )
                h.setpointer(pointer, 'pka', seg.cal12 )
                h.setpointer(pointer, 'pka', seg.cal13 )
                h.setpointer(pointer, 'pka', seg.kir )
                
                
                if sec.name().find('soma') >= 0:
                    
                    
                    # can is only distributed to the soma section
                    h.setpointer(pointer, 'pka', seg.can )
                    
                    
                    
                    
    # synaptic modulation ================================================================
    if sim == 'synMod':
        
        # draw random modulation factors for synapses from
        # intervals given by range_list[[min,max]], where 1 is no modulation.  
        glut_f, glut_f_norm     = set_rand_synapse(['amp', 'nmd'], base_mod, max_mod,   \
                                                    range_list=[[0.9,1.6], [0.9,1.6]]   )
                                                    
        gaba_f, gaba_f_norm     = set_rand_synapse(['gab'],        base_mod, max_mod,   \
                                                    range_list=[[0.6,1.4]]              )
        
        syn_fact = glut_f + gaba_f
        
            
        I_d = {}
        ns  = {}
        nc  = {}
        Syn = {}
        
        for sec in h.allsec():
            
            # create one glutamatergic and one gabaergic synapse per section
            if sec.name().find('dend') >= 0:
                
                
                # create a glut synapse
                make_random_synapse(ns, nc, Syn, sec, 0.5,          \
                                        NS_interval=1000.0/17.0,    \
                                        NC_conductance=0.15e-3,     \
                                        S_tau_dep=100               )
                                        
                                        
                # create a gaba synapse
                make_random_synapse(ns, nc, Syn, sec, 0.0,          \
                                        Type='tmgabaa',             \
                                        NS_interval=1000.0/4.0,     \
                                        NC_conductance=0.45e-3      )
                
                
                # set pointer(s)
                h.setpointer(pointer, 'pka', Syn[sec])
                h.setpointer(pointer, 'pka', Syn[sec.name()+'_gaba'])
                
                
                # configure
                Syn[sec].base       = base_mod
                Syn[sec].f_ampa     = glut_f_norm[0]
                Syn[sec].f_nmda     = glut_f_norm[1]
                
                Syn[sec.name()+'_gaba'].base    = base_mod
                Syn[sec.name()+'_gaba'].f_gaba  = gaba_f_norm[0]
            
            
            elif sec.name().find('axon') >= 0: 
                
                # don't modulate segments in the axon
                continue 
                  
            
            # set modulation of channels
            for seg in sec:
                
                for mech in seg:
                    
                    
                    # turn of or set modulation
                    if mech.name() in not2mod:
                        
                        mech.factor = 0.0
                        print(mech.name(), 'and channel:', not2mod, mech.factor, sec.name())
                        
                    elif mech.name() in mod_list:
                    
                        mech.base       = base_mod
                        index           = mod_list.index( mech.name() )
                        mech.factor     = factors[index]
                        
                    
                    
                    
                    
    
    # solver------------------------------------------------------------------------------            
    cvode = h.CVode()
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop:
    
        if dynMod == 1:
        
            if h.t > da_tstart: 
                
                casc.DA = alpha(da_tstart, da_peak, da_tau) 
                
        h.fadvance()
        
    
    
    # save output ------------------------------------------------------------------------
    
    ID          = ''
    all_factors = mod_fact + syn_fact
    
    for i,mech in enumerate(mod_list+['amp', 'nmd', 'gab']):
        ID = ID + mech + str( int(all_factors[i]*100) )
    
    if dynMod == 1:
        # DA transient
        save_vector(tm, vm, ''.join(['./Results/Dynamic/spiking_', str(run), '_', ID, '.out']) )
        
    else: 
        # no transient
        save_vector(tm, vm, ''.join(['./Results/Dynamic/spiking_', str(run), '_control.out']) )
        
    
    if run == 0:
        # save substrate concentrations
        names = ['Target1p', 'cAMP', 'Gbgolf', 'D1RDAGolf']
        
        for i,substrate in enumerate([pka, camp, gbg, gprot]):
            save_vector(tm, substrate, './Results/Dynamic/substrate_'+names[i]+'.out' )
        
    
    
    
    # save current from the synapses
    #for key in I_d:
        #save_vector(tm, I_d[key], ''.join(['./I_', key, '_', str(run), '.out']) )
            
    # when run large scale (on supercomputer) all iterarations from one core 
    # were stored in one dict to decrease the number of files. Only factors and spikes
    # were recorded        
    #spikes      = fun.getSpikedata_x_y(tm,vm) 
    #RES[run]    = {'factors': mod_fact + syn_fact, 'spikes': spikes}
        
                



# Start the simulation.
# Function needed for HBP compability  ===================================================
if __name__ == "__main__":
    
    
    print('starting sim')
    
    
    n_runs_control      = 5
    n_runs_modulated    = 5
    
    # other targets: ['Target1p', 'cAMP', 'Gbgolf',  'D1RDAGolf']
    
    # modulated
    for n in range(n_runs_modulated):
        main( par="./params-rob.json",          \
                    amp=0.0,                    \
                    run=n,                      \
                    simDur=2000,                \
                    stimDur=3000,               \
                    sim='synMod',               \
                    dynMod=1,                   \
                    target='Target1p',          \
                    not2mod=[]                  )
    
    # control (dynMod = 0)
    for n in range(n_runs_modulated,n_runs_modulated+n_runs_modulated):
        main( par="./params-rob.json",          \
                    amp=0.0,                    \
                    run=n,                      \
                    simDur=2000,                \
                    stimDur=3000,               \
                    sim='synMod',               \
                    dynMod=0,                   \
                    target='Target1p',          \
                    not2mod=[]                  )
        
        '''
        # write result to file every 20 iterations (if run large scale in parallel)
        if n % 20 == 0:  
            
            print('in save loop') 
            
            mod_list    = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can']
            syn_mod     = ['amp', 'nmd', 'gab'] 
            mod_fact    = RES[0]['factors']
            ID          = ''
            
            for i,mech in enumerate(mod_list+syn_mod):
                ID = ID + mech + str( int(mod_fact[i]*100) )
            
            fun.save_obj(RES, ''.join(['ModAll_', target, '-', ID]) )
        '''
            
            
            
    print('plotting')            
    fun.plot_fig6B('./Results/Dynamic/')            
    
                        
                                                    
    
                                                    
                                                    
                                                    
                                                    
    
    
    
          
    
        

