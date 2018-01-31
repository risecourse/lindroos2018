#
# MSN model
#
# Robert Lindroos
# 
# Original version by 
# Alexander Kozlov <akozlov@kth.se>
# Kai Du <kai.du@ki.se>
#


from __future__ import print_function, division
import mpi4py as MPI
import sys
import json
import numpy as np
from neuron import h
from math import exp
import pickle
import os



mod        = "./mod/"
defparams  = "./params-rob.json"
morphology = "./morphology/"

with open('./substrates.json') as file:
    SUBSTRATES = json.load(file)



#h.nrn_load_dll('/pdc/vol/neuron/7.4-py27/x86_64/.libs/libnrnmech.so')
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')



# global result dict
RES = {}


# Distributions:
'''
T-type Ca: (1.0/(1+np.exp((v-70)/-4.5)))
naf (den): (0.1 + 0.9/(1 + np.exp((v-60.0)/10.0)))

'''

#                        1,      0.8, 2.5,60,-25
#calculate_distribution(d3, dist, a4, a5, a6, a7, g8)
def calculate_distribution(d3, dist, a4, a5, a6, a7, g8):
    # d3 is the distribution type:
    #     0 linear, 1 sigmoid, 2 exponential
    #     3 step for absolute distance (in microns)
    # dist is the somatic distance
    # a4-7 is distribution parameters 
    # g8 is the maximal conductance
    if   d3 == 0: 
        value = a4 + a5*dist
    elif d3 == 1: 
        value = a4 + a5/(1 + exp((dist-a6)/a7) )
    elif d3 == 2: 
        value = a4 + a5*exp((dist-a6)/a7)
    elif d3 == 3:
        if (dist > a6) and (dist < a7):
            value = a4
        else:
            value = a5
            
    if value < 0:
        value = 0
        
    value = value*g8
    return value 


  
    
def alpha(tstart, gmax, tau):
    # calc and returns a "magnitude" using an alpha function -> used for [DA] in cascade
    
    v = 1 - (h.t - tstart) / tau
    e = exp(v)
    mag = gmax * (h.t - tstart) / tau * e
    
    return mag
    
    
    
def save_vector(t, v, outfile):
    
    with open(outfile, "w") as out:
        for time, y in zip(t, v):
            out.write("%g %g\n" % (time, y))
            
            
            
def calc_rand_Modulation(mod_list, range_list=False, distribution='centered'):
    
    '''
    calc random modulation values between 0 and 2 (i.e. max +/- 100%).
    
    Values close to 1 (no mod) are given higer probability.
    this is achived by drawing two uniform random numbers between 0 and 1, and subtracts one from the other.
    This gives values ranging from -1 to 1. By adding one we end upp at the wanted value.
    
    If a range_list is supplied the range of the values can be shifted from [0,2] to any
    other range. The range list must have the same length as mod_list and hold lists of
    [min, max] values.
    
    '''
    
    mod_factors = []
    
    A=0
    B=2
    
    for i,channel in enumerate(mod_list):
        
        if distribution=='centered':
            factor = 1.0 + ( np.random.uniform() - np.random.uniform() )
        elif distribution=='inv_centered':
            factor = 1.0 + ( np.random.uniform() - np.random.uniform() )
            if factor <= 1:
                factor = factor + 1.0
            else:
                factor = factor - 1.0
        elif distribution=='uniform':
            factor = 2.0 * np.random.uniform()
        else:
            print('Error in MF distribution--line ~121')
            eegsjsd
        
        
        if range_list:
            
            a       = range_list[i][0]
            b       = range_list[i][1]
            
            factor = (b-a) / (B-A) * (factor-A) + a
       
        mod_factors.append(factor)
        
    return mod_factors 
    
    

# 'save/'+ 
def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open('save/' + name + '.pkl', 'rb') as f:
        return pickle.load(f) 
        
        
        
def getSpikedata_x_y(x,y):
    
    ''' 
    There's probably a Neuron function for this--use instead?
    
    getSpikedata_x_y(x,y) -> return count
    
    Extracts and returns the number of spikes from spike trace data.
    
    # arguments
    x = time vector
    y = vm vector 
    
    # returns
    count = number of spikes in trace (int)
    
    # extraction algorithm
    -threshold y and store list containing index for all points larger than 0 V 
    -sorts out and counts the index that are the first one(s) crossing the threshold, i.e. 
        the first index of each spike. This is done by looping over all index and check if 
        the index is equal to the previous index + 1. If not it is the first index of a 
        spike.
        
        If no point is above threshold in the trace the function returns 0.
        
    '''
    
    count = 0
    spikes = []
    
    # pick out index for all points above zero potential for potential trace
    spikeData = [i for i,v in enumerate(y) if v > 0]

    # if no point above 0
    if len(spikeData) == 0:
        
        return spikes
	
    else:
        # pick first point of each individual transient (spike)...
        for j in range(0, len(spikeData)-1):
            if j==0:
                
                count += 1
                spikes.append(x[spikeData[j]])
		
            # ...by checking above stated criteria
            elif not spikeData[j] == spikeData[j-1]+1:
                count += 1
                spikes.append(x[spikeData[j]])
            
    return spikes 


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
    
    syn_fact = calc_rand_Modulation(channel_list, range_list=range_list, distribution='uniform')
        
    # normalize factors to max-value of pointer substrate
    normalized_factors     = []
    for i,mech in enumerate(channel_list):
        
        normalized_factors.append( (syn_fact[i] - 1) / (max_mod - base_mod)  ) 
        
    return syn_fact, normalized_factors     
 
# ======================= the MSN class ==================================================

class MSN:
    def __init__(self, params=defparams, factors=None):
        Import = h.Import3d_SWC_read()
        Import.input(morphology + 'latest_WT-P270-20-14ak.swc')
        imprt = h.Import3d_GUI(Import, 0)
        imprt.instantiate(None)
        h.define_shape()
        # h.cao0_ca_ion = 2  # default in nrn
        h.celsius = 35
        self._create_sectionlists()
        self._set_nsegs()
        self.v_init = -80
        for sec in self.allseclist:
            sec.Ra = 150
            sec.cm = 1.0
            sec.insert('pas')
            #sec.g_pas = 1e-5 # set using json file
            sec.e_pas = -70 # -73
        for sec in self.somalist:
            sec.insert('naf')
            sec.insert('kaf')
            sec.insert('kas')
            sec.insert('kdr')
            sec.insert('kir')
            sec.ena = 50
            sec.ek = -85 # -90
            sec.insert('cal12')
            sec.insert('cal13')
            sec.insert('car')
            sec.insert('cadyn')
            sec.insert('caldyn')
            sec.insert('sk')
            sec.insert('bk')
            sec.insert('can')
	    #sec.kb_cadyn = 200.
        for sec in self.axonlist:
            sec.insert('naf')
            #sec.insert('kaf')
            sec.insert('kas')
            #sec.insert('kdr')
            #sec.insert('kir')
            sec.ena = 50
            sec.ek = -85 # -90
        for sec in self.dendlist:
            sec.insert('naf')
            sec.insert('kaf')
            sec.insert('kas')
            sec.insert('kdr')
            sec.insert('kir')
            sec.ena = 50
            sec.ek = -85 # -90
            sec.insert('cal12')
            sec.insert('cal13')
            sec.insert('car')
            sec.insert('cadyn')
            sec.insert('caldyn')
            sec.insert('sk')
            sec.insert('bk')
            sec.insert('cat32')
            sec.insert('cat33')

        with open(params) as file:
            par = json.load(file)

        self.distribute_channels("soma", "g_pas", 0, 1, 0, 0, 0, float(par['g_pas_all']['Value']))
        self.distribute_channels("axon", "g_pas", 0, 1, 0, 0, 0, float(par['g_pas_all']['Value']))
        self.distribute_channels("dend", "g_pas", 0, 1, 0, 0, 0, float(par['g_pas_all']['Value']))

        self.distribute_channels("soma", "gbar_naf", 0, 1, 0, 0, 0, float(par['gbar_naf_somatic']['Value']),factors=factors)
        self.distribute_channels("soma", "gbar_kaf", 0, 1, 0, 0, 0, float(par['gbar_kaf_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kas", 0, 1, 0, 0, 0, float(par['gbar_kas_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kdr", 0, 1, 0, 0, 0, float(par['gbar_kdr_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kir", 0, 1, 0, 0, 0, float(par['gbar_kir_somatic']['Value']))
        self.distribute_channels("soma", "gbar_sk", 0, 1, 0, 0, 0, float(par['gbar_sk_somatic']['Value']))
        self.distribute_channels("soma", "gbar_bk", 0, 1, 0, 0, 0, float(par['gbar_bk_somatic']['Value']))
        
        #self.distribute_channels("axon", "gbar_naf", 0, 1, 0, 0, 0, float(par['gbar_naf_somatic']['Value']),factors=factors)
        self.distribute_channels("axon", "gbar_naf", 3, 1, 1.1, 30, 500, float(par['gbar_naf_axonal']['Value']),factors=factors)
        #self.distribute_channels("axon", "gbar_naf", 3, 1, 1.1, 20, 500, float(par['gbar_naf_axonal']['Value']))
        #self.distribute_channels("dend", "gbar_naf", 1, 1,  1.2, 30, -5, float(par['gbar_naf_axonal']['Value']))
        self.distribute_channels("axon", "gbar_kas", 0, 1, 0, 0, 0, float(par['gbar_kas_axonal']['Value']))
        
        self.distribute_channels("dend", "gbar_naf", 1, 0.1, 0.9, 60.0, 10.0, float(par['gbar_naf_basal']['Value']),factors=factors)
        self.distribute_channels("dend", "gbar_kaf", 1, 1,  0.5, 120, -30, float(par['gbar_kaf_basal']['Value']))
        #self.distribute_channels("dend", "gbar_naf", 0, 1, -0.0072, 0, 0, float(par['gbar_naf_basal']['Value']))
        #self.distribute_channels("dend", "gbar_kaf", 0, 1,  0.0167, 0, 0, float(par['gbar_kaf_basal']['Value']))
        self.distribute_channels("dend", "gbar_kas", 2, 1, 9.0, 0.0, -5.0, float(par['gbar_kas_basal']['Value']))
        self.distribute_channels("dend", "gbar_kdr", 0, 1, 0, 0, 0, float(par['gbar_kdr_basal']['Value']))
        self.distribute_channels("dend", "gbar_kir", 0, 1, 0, 0, 0, float(par['gbar_kir_basal']['Value']))
        self.distribute_channels("dend", "gbar_sk", 0, 1, 0, 0, 0, float(par['gbar_sk_basal']['Value']))
        self.distribute_channels("dend", "gbar_bk", 0, 1, 0, 0, 0, float(par['gbar_bk_basal']['Value']))

        self.distribute_channels("soma", "pbar_cal12", 0, 1, 0, 0, 0, 1e-5)
        self.distribute_channels("soma", "pbar_cal13", 0, 1, 0, 0, 0, 1e-6)
        self.distribute_channels("soma", "pbar_car", 0, 1, 0, 0, 0, 1e-4)
        self.distribute_channels("soma", "pbar_can", 0, 1, 0, 0, 0, 3e-5)
        #self.distribute_channels("soma", "kb_cadyn", 0, 1, 0, 0, 0, 200.0)
        self.distribute_channels("dend", "pbar_cal12", 0, 1, 0, 0, 0, 1e-5)
        self.distribute_channels("dend", "pbar_cal13", 0, 1, 0, 0, 0, 1e-6)
        self.distribute_channels("dend", "pbar_car", 0, 1, 0, 0, 0, 1e-4)
        self.distribute_channels("dend", "pbar_cat32", 1, 0, 1.0, 70.0, -4.5, 1e-7)
        self.distribute_channels("dend", "pbar_cat33", 1, 0, 1.0, 70.0, -4.5, 1e-8)

    def _create_sectionlists(self):
        self.allsecnames = []
        self.allseclist = h.SectionList()
        for sec in h.allsec():
            self.allsecnames.append(sec.name())
            self.allseclist.append(sec=sec)
        self.nsomasec = 0
        self.somalist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('soma') >= 0:
                self.somalist.append(sec=sec)
                if self.nsomasec == 0:
                    self.soma = sec
                self.nsomasec += 1
        self.axonlist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('axon') >= 0:
                self.axonlist.append(sec=sec)
        self.dendlist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('dend') >= 0:
                self.dendlist.append(sec=sec)

    def _set_nsegs(self):
        for sec in self.allseclist:
            sec.nseg = 2*int(sec.L/40.0)+1
        for sec in self.axonlist:
            sec.nseg = 2  # two segments in axon initial segment

    def _max_dist(self, axon_excluding=True):
        h.distance(sec=self.soma)
        dmax = 0
        for sec in self.allseclist:
	        if axon_excluding and sec.name().find('axon') == 0: continue
                dmax = max(dmax, h.distance(1, sec=sec))
        return dmax

    def distribute_channels(self, as1, as2, d3, a4, a5, a6, a7, g8, factors=None):
        h.distance(sec=self.soma)
        dmax = self._max_dist()
        for sec in self.allseclist:
            if sec.name().find(as1) >= 0:
                for seg in sec:
                    dist = h.distance(seg.x, sec=sec)
                    val = calculate_distribution(d3, dist, a4, a5, a6, a7, g8)
                    cmd = 'seg.%s = %g' % (as2, val)
                    exec(cmd)
                    
                    #names = ['mVhalf_naf', 'hVhalf_naf', 'mSlope_naf', 'hSlope_naf']
                    names  = ['taum_naf', 'tauh_naf', 'taun_naf']
                    
                    if factors:
                        for i,factor in enumerate(factors):
                            
                            print('updating factors! ', factors, factor, names[i])
                            
                            # add minus here if running with slope and half activation; cmd = 'seg.%s = %g' % (names[i], -factor)
                            cmd = 'seg.%s = %g' % (names[i], factor)
                            exec(cmd)
                            
 
                  
   
#=========================================================================================



def main(par="./params-msn.json", \
                            sim='vm',       \
                            amp=0.265,      \
                            run=None,       \
                            modulation=1,   \
                            simDur=7000,    \
                            stimDur=900,    \
                            factors=None,   \
                            section=None,   \
                            randMod=None,   \
                            testMode=False, \
                            target=None,    \
                            chan2mod=['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can'] ): 
    
    
    
    print(locals())
    
    # initiate cell
    cell = MSN(params=par, factors=factors)
        
    # set cascade ---- move to MSN def?
    casc = h.D1_reduced_cascade2_0(0.5, sec=cell.soma) # other cascades also possible...
    
    if target:
        cmd = 'pointer = casc._ref_'+target
        exec(cmd)
        
        base_mod    = SUBSTRATES[target][0]
        max_mod     = SUBSTRATES[target][1]
        
    else:
        pointer     = casc._ref_Target1p    #Target1p   #totalActivePKA    (if full cascade used)
        base_mod    = casc.init_Target1p
        max_mod     = 2317.1
    
    # cAMP; init: 38.186016
    
    # set edge of soma as reference for distance 
    h.distance(1, sec=h.soma[0])
    
    # set current injection
    stim = h.IClamp(0.5, sec=cell.soma)
    stim.amp = amp  
    stim.delay = 100
    stim.dur = stimDur            # 2ms 2nA to elicit single AP, following Day et al 2008 Ca dyn    
    
    # record vectors
    tm = h.Vector()
    tm.record(h._ref_t)
    vm = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    # substrates
    pka = h.Vector()
    pka.record(casc._ref_Target1p)
    camp = h.Vector()
    camp.record(casc._ref_cAMP)
    gprot = h.Vector()
    gprot.record(casc._ref_D1RDAGolf) #D1RDAGolf
    gbg   = h.Vector()
    gbg.record(casc._ref_Gbgolf) #Gbgolf
    
    # peak n dipp parameters
    da_peak   = 500   # concentration [nM]
    da_tstart = 1000    # stimulation time [ms]
    da_tau    = 500    # time constant [ms]
    
    
    tstop = simDur               # [ms]
    
    
    # all channels to modulate
    mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can' ]
    
    
    not2mod = [] #['kaf']
    
    
    # find channels that should not be modulated
    for chan in mod_list:
        
        if chan not in chan2mod:
            
            not2mod.append(chan)
    
    
    # for random modulation: modValues = np.arange(0.1, 2.0, 0.1) -------------------------
    if randMod == 1:
        
        # new factors every run
        mod_fact = calc_rand_Modulation(mod_list, range_list=[[0.60,0.80],    \
                                                              [0.65,0.85],  \
                                                              [0.75,0.85],  \
                                                              [0.85,1.25],  \
                                                              [1.0,2.0],    \
                                                              [1.0,2.0],    \
                                                              [0.0,1.0]],
                                                              distribution='uniform'  )
        
        # keep old factors
        '''
        if run == 0:
        
            mod_fact = calc_rand_Modulation(mod_list)
            
        else:
            
            mod_fact = RES['factors']
        '''
    
    else:
        mod_fact = [ 0.8, 0.8, 0.8, 1.25, 2.0, 2.0, 0.5  ]
           
    
    # noormalize factors to  target values seen in simulation
    factors     = []
    for i,mech in enumerate(mod_list):
        
        factor  = (mod_fact[i] - 1) / (max_mod - base_mod) #2317.1
        
        factors.append(factor)
        
        #print(mech, mod_fact[i], factor) # --------------------------------------------------------
            
    
    # set pointers 
    for sec in h.allsec():
        
        for seg in sec:
            
            # naf and kas is in all sections
            h.setpointer(pointer, 'pka', seg.kas )
            h.setpointer(pointer, 'pka', seg.naf )
            
            if sec.name().find('axon') < 0:    
                
                # these channels are not in the axon section
                
                h.setpointer(pointer, 'pka', seg.kaf )
                h.setpointer(pointer, 'pka', seg.cal12 )
                h.setpointer(pointer, 'pka', seg.cal13 )
                h.setpointer(pointer, 'pka', seg.kir )
                #h.setpointer(pointerc, 'pka', seg.car )
                
                if sec.name().find('soma') >= 0:
                    
                    # can is only distributed to the soma section
                    h.setpointer(pointer, 'pka', seg.can )
                    
                    
                    
                    


    # synaptic modulation ================================================================
    if sim == 'synMod':
        
        
        # draw random modulation factors (intervals given by range_list[[min,max]]  
        glut_f, glut_f_norm     = set_rand_synapse(['amp', 'nmd'], base_mod, max_mod,   \
                                                    range_list=[[0.9,1.6], [0.9,1.6]]   )
                                                    
        gaba_f, gaba_f_norm     = set_rand_synapse(['gab'],        base_mod, max_mod,   \
                                                    range_list=[[0.6,1.4]]              )
        
        syn_fact = glut_f + gaba_f
            
        I_d={}
        
        ns = {}
        nc = {}
        Syn = {}
        for sec in h.allsec():
            if sec.name().find('dend') >= 0:
                
                # create a glut synapse
                make_random_synapse(ns, nc, Syn, sec, 0.5,          \
                                        NS_interval=1000.0/20.0,    \
                                        NC_conductance=0.165e-3,     \
                                        S_tau_dep=100               )
                                        
                # create a gaba synapse
                make_random_synapse(ns, nc, Syn, sec, 0.0,          \
                                        Type='tmgabaa',             \
                                        NS_interval=1000.0/5.0,     \
                                        NC_conductance=0.495e-3      )
                
                # set pointer(s)
                h.setpointer(pointer, 'pka', Syn[sec])
                h.setpointer(pointer, 'pka', Syn[sec.name()+'_gaba'])
                
                # set (random?) modulation
                Syn[sec].base    = base_mod
                
                #randMod?
                if randMod == 1:
                    Syn[sec].f_ampa     = glut_f_norm[0]
                    Syn[sec].f_nmda     = glut_f_norm[1]
                else:
                    Syn[sec].f_ampa     = 0
                    Syn[sec].f_nmda     = 0
                
                if randMod == 1:
                    Syn[sec.name()+'_gaba'].base    = base_mod
                    Syn[sec.name()+'_gaba'].f_gaba  = gaba_f_norm[0]
                else:
                    Syn[sec.name()+'_gaba'].f_gaba  = 0
                    
                '''
                # record synaptic current from synapse

                I_d[sec.name()]             = h.Vector()

                I_d[sec.name()].record(Syn[sec]._ref_i)
                
                I_d[sec.name()+'_gaba']     = h.Vector()
                I_d[sec.name()+'_gaba'].record(Syn[sec.name()+'_gaba']._ref_i)
                '''
            
            
            elif sec.name().find('axon') >= 0: 
                continue   
            
            if randMod == 1:
                for seg in sec:
                    
                    for mech in seg:
                        
                        if mech.name() in not2mod:
                            
                            mech.factor = 0.0
                            print(mech.name(), 'and channel:', not2mod, mech.factor, sec.name())
                            
                        elif mech.name() in mod_list:
                        
                            mech.base       = base_mod
                            index           = mod_list.index( mech.name() )
                            mech.factor     = factors[index]
                    
                    
                    
                    
    
    # dynamical modulation        
    elif sim == 'modulation':
        
        print('inne ', sim)
        
        for sec in h.allsec():
            
            for seg in sec:
                
                for mech in seg:
                    
                    # if this first statement is active the axon will not be modulated
                    '''if sec.name().find('axon') >= 0     \
                            and mech.name() in mod_list:
                            
                        mech.factor = 0.0
                        print(sec.name(), seg.x, mech.name() )'''
                        
                    if mech.name() in not2mod:
                        
                        mech.factor = 0.0
                        print(mech.name(), 'and channel:', not2mod, mech.factor, sec.name())
                        
                    elif mech.name() in mod_list:
                    
                        mech.base       = base_mod
                        index           = mod_list.index( mech.name() )
                        mech.factor     = factors[index]
                    
                    
                        
                    
    
    # static modulation
    elif sim == 'directMod':
        
        print('inne ', sim)
    
        for sec in h.allsec():
            
            for seg in sec:
                
                for mech in seg:
                
                    if mech.name() in mod_list: 
                        
                        if sec.name().find('axon') < 10: # 0 no axon modulated; 10 all sections
                        
                            factor = mod_fact[mod_list.index(mech.name() )]
                            
                            if mech.name() in not2mod:
                                mech.factor = 0.0
                            elif mech.name()[0] == 'c':
                                pbar = mech.pbar
                                mech.pbar = pbar * factor
                                #print(''.join(['setting pbar ', mech.name(), ' ', str(factor) ]) )
                            else:
                                gbar = mech.gbar
                                mech.gbar = gbar * factor
                                #if seg.x < 0.2:
                                #print(''.join(['setting gbar ', mech.name(), ' ', str(factor), ' ', str(gbar), ' ', str(factor*gbar) ]) )
                        
                        else:
                            
                            print(sec.name(), seg.x, sec.name().find('axon'))
    
    
    
    # solver------------------------------------------------------------------------------            
    cvode = h.CVode()
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop:
    
        if modulation == 1:
        
            if h.t > da_tstart: 
                
                # set DA and ACh values (using alpha function)
                casc.DA = alpha(da_tstart, da_peak, da_tau) 
                #casc.ACh = ach_base - alpha(ach_tstart, ach_base, ach_tau)
                
        h.fadvance()
        
    
    
    # save output
    if sim in ['vm', 'directMod', 'modulation', 'synMod']:
        
        if testMode:
        
            mod_fact = mod_fact + syn_fact
                
            ID = ''
            
            for i,mech in enumerate(mod_list+syn_mod):
                
                ID = ID + mech + str( int(mod_fact[i]*100) )
                            
            save_vector(tm, vm, ''.join(['./spiking_', str(run), '_', ID, '.out']) )
            '''
            names = ['Target1p', 'cAMP', 'Gbgolf', 'D1RDAGolf']
            data  = {}
            for i,substrate in enumerate([pka, camp, gbg, gprot]):
                save_vector(tm, substrate, './substrate_'+names[i]+'.out' )
                print(min(substrate), max(substrate))
                data[names[i]] = [min(substrate), max(substrate)] 
            
            print(ID)
            with open('substrates.json', 'w') as outfile:
                json.dump(data, outfile)'''
            
            #for key in I_d:
                
                #save_vector(tm, I_d[key], ''.join(['./I_', key, '_', str(run), '.out']) )
                
                
        spikes      = getSpikedata_x_y(tm,vm) 
        
        RES[run]    = {'factors': mod_fact + syn_fact, 'spikes': spikes}
        
                



# if run from terminal...   ===============================================================
if __name__ == "__main__":
    
    sys.argv = ['a', 'synMod']
    
    print('starting sim')
    
    factors = [1] #list(itertools.product(taum, tauh, taun))
    
    if sys.argv[1] in ['vm', 'directMod']:  #for current in currents:
    
        mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can']
        
        randMod = 1
        
        currents = [0] #np.arange(300,360,10)
        
        for n in range(200):
        
            RES = {}
        
            for i,current in enumerate(currents):
                main( par="./params-rob.json",      \
                        amp=current*1e-3,           \
                        run=i,                      \
                        modulation=0,               \
                        simDur=3000,                \
                        stimDur=3000,               \
                        sim=sys.argv[1],            \
                        randMod=randMod,            \
                        chan2mod=mod_list           )
            
            
            mod_fact = RES['factors']
            
            ID = ''
            
            for i,mech in enumerate(mod_list):
                
                ID = ID + mech + str( int(mod_fact[i]*100) )
        
            save_obj(RES, ''.join(['DRM-', ID]) )
                
                
                
    elif sys.argv[1] in 'synMod':  #for current in currents:
    
        mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can']
        syn_mod  = ['amp', 'nmd', 'gab']
        
        currents = [0] #np.arange(300,360,10)
        
        RES = {}
        
        testMode = False
        
        if testMode:
            n_runs = 4
        else:
            n_runs = 401
            
        # ['Target1p', 'cAMP', 'Gbgolf',  'D1RDAGolf']
        for n in range(n_runs):
            main( par="./params-rob.json",          \
                        amp=0.0,                    \
                        run=n,                      \
                        simDur=2000,                \
                        stimDur=3000,               \
                        sim=sys.argv[1],            \
                        modulation=1,               \
                        randMod=1,                  \
                        testMode=testMode,          \
                        target='Target1p',          \
                        chan2mod=mod_list           )
            
            
            if not testMode and n % 20 == 0:  
                
                print('in save loop')  
                            
                mod_fact = RES[0]['factors']
                
                ID = ''
                
                for i,mech in enumerate(mod_list+syn_mod):
                    
                    ID = ID + mech + str( int(mod_fact[i]*100) )
                
                save_obj(RES, ''.join(['HighExcUniModAll_target1-', ID]) )
                        
                                                    
    
                                                    
                                                    
                                                    
                                                    
    
    
    
          
    
        

