#
'''
MSN model used in Lindroos et al., (2018). Frontiers

Robert Lindroos (RL) <robert.lindroos at ki.se>
 
The MSN class and most channels were implemented by 
Alexander Kozlov <akozlov at kth.se>
with updates by RL

Implemented in colaboration with Kai Du <kai.du at ki.se>
'''




from neuron import h
import numpy                as np
import plot_functions       as fun
import MSN_builder          as build
import os.path


h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')



# global result dict
RES = {}


   
    
def save_vector(t, v, outfile):
    
    with open(outfile, "w") as out:
        for time, y in zip(t, v):
            out.write("%g %g\n" % (time, y))
            
    
            
            
def calc_rand_Modulation(mod_list, range_list=False):
    '''
    uses numpy to draw random modulation factors in range [0,2],
    from a uniform distribution, for each channel in mod_list.
    
    The factors can also be linearly mapped to an arbitrary interval. 
    This is done if a range_list is given.
    
    mod_list     = list of channels to be modulated
    range_list   = list of [min, max] values to be used in modulation. 
                    Must have same length as mod_list.
    '''
    # TODO: Figure out how treatments would affect the dopamine concentration
    # Try doubling the modulation factor MF (or some other increase) to represent
    # the effect of the therapeutics.
    #
    # 
    
    
    MF = 2.0 # Originally 2.0, we could alter this to represent increased
             # or decreased dopamine. We could also alter the range of
             # modulation values that can be present in the model, perhaps
             # the neurons are more consistently modulated in the presence
             # of some of these therapeutics. Ex:
                 # MAO-B Inhibitor:  A =    B =      MF =
                 # Dopamine Agonist: A =    B =      MF =
                 # Levodopa:         A =    B =      MF =
             
    mod_factors = []
    
    A=0
    B=2
    
    for i,channel in enumerate(mod_list):
        
        factor = MF * np.random.uniform()
    
        if range_list:
            
            a       = range_list[i][0]
            b       = range_list[i][1]
            
            factor = (b-a) / (B-A) * (factor-A) + a
       
        mod_factors.append(factor)
        
    return mod_factors 
        
         



   
#=========================================================================================


# in the dynamimcal modulation, the channels are connected to one substrate of the cascade.
# base modulation (control) is assumed for base value of the substrate and full modulation
# is assumed when the substrate level reaches its maximal value. Linear scaling is used 
# between these points.
def main(par="./params-msn.json", \
                            sim='vm',       \
                            amp=0.265,      \
                            run=None,       \
                            simDur=7000,    \
                            stimDur=900,    \
                            not2mod=[],     \
                            modulate_axon=False): 
    
    
    print('-iter:', run, '-amp:', amp, '[nA]')
    
    
    
    # initiate cell
    cell    =   build.MSN(  params=par,                             \
                            morphology='latest_WT-P270-20-14ak.swc' )
    
    
    # set cascade--not connected to channels in this script, 
    # but used for setting pointers needed in the channel mechnisms
    casc    =   h.D1_reduced_cascade2_0(0.5, sec=cell.soma) 
    pointer =   casc._ref_Target1p    
       
    
    # set edge of soma as reference for dendritic distance 
    h.distance(1, sec=h.soma[0])
    
    
    # set current injection
    stim        =   h.IClamp(0.5, sec=cell.soma)
    stim.amp    =   amp  
    stim.delay  =   100
    stim.dur    =   stimDur   
    
    
    # record vectors
    tm = h.Vector()
    tm.record(h._ref_t)
    vm = h.Vector()
    vm.record(cell.soma(0.5)._ref_v) 
     
    
    tstop       = simDur 
    # dt = default value; 0.025 ms (25 us)
    
    
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
                    
                    
                    
                    
    # static modulation ================================================================
    if sim == 'directMod':
        
        
        # channels to modulate channels
        mod_list    = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can' ] 
          
        
        if amp == 0.32:
            
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
                                                                  
            print('factors drawn:', mod_fact)
            
        else:
            
            # use same factors for all three current injections
            mod_fact = RES[run]['factors']
            
            
    
        for sec in h.allsec():
            
            # modulate axon?
            if sec.name().find('axon') >= 0:
                if not modulate_axon:
                    continue
            
            for seg in sec:
                
                for mech in seg:
                    
                    if mech.name() in not2mod:
                        continue
                    
                    elif mech.name() in mod_list: 
                        
                        # get factor from list
                        factor = mod_fact[mod_list.index(mech.name() )]
                        
                        if mech.name()[0] == 'c':
                            # Ca channels
                            pbar        = mech.pbar
                            mech.pbar   = pbar * factor
                            
                        else:
                            # non Ca channels
                            gbar        = mech.gbar
                            mech.gbar   = gbar * factor
                        
                            
                        
                    
                    
    # solver------------------------------------------------------------------------------
                
    cvode = h.CVode()
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop: 
                
        h.fadvance()
        
    
    # save output ------------------------------------------------------------------------
    # Extract spikes and save together with modulation factors
    
    spikes      = fun.getSpikedata_x_y(tm,vm) 
    amplitude   = int(amp * 1000)
    
    
    if amp == 0.32:
    
        if run == 0:
        
            RES[run]    = {amplitude: spikes}
            
        else:
        
            RES[run]    = {'factors': mod_fact, amplitude: spikes}
            
    else:
        
        RES[run][amplitude] = spikes  
        
                



# Start the simulation.
# Function needed for HBP compability  ===================================================
if __name__ == "__main__":
    
    
    simulate    = True
    file_path   = 'Results/static_spikes.pkl'
    
    
    if simulate:
        
        
        print('starting sim')
        
        mod_list    = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can']
        currents    = np.arange(320,345,10)
        n_start     = 1
        n_runs      = 10
        
        
        # check if previous simulation exists
        if os.path.exists(file_path):
            
            # load file
            RES = fun.load_obj( 'Results/static_spikes.pkl' )
            
            # find max of key in RES
            for key in RES:
                
                if key > n_start:
                    n_start = key
            
            # add one to the maximum (not to overwrite one trace in the result dict)
            n_start += 1
        
        
        # control (without modulation; only if not run before)
        if n_start == 1:
            for current in currents:
                main( par="./params_dMSN.json",      \
                        amp=current*1e-3,           \
                        run=0,                      \
                        simDur=1000,                \
                        stimDur=900,                \
                        sim='vm',                   \
                        not2mod=[],                 \
                        modulate_axon=False         )
        
        
        # randomly modulated
        for n in range(n_start,n_start+n_runs):
        
            print('\n----------------------')
            
            for current in currents:
                main( par="./params_dMSN.json",      \
                        amp=current*1e-3,           \
                        run=n,                      \
                        simDur=1000,                \
                        stimDur=900,                \
                        sim='directMod',            \
                        not2mod=[],                 \
                        modulate_axon=False         )
            
        fun.save_obj( RES, file_path.split('.')[0] )
        
    else:
    
        RES = fun.load_obj( file_path )       
    
        
    print('plotting')    
    fun.plot_fig4C(RES)           
    
                        
                    
                    
                        
                    
    
    
