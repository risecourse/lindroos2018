#
# MSN model
#
# Robert Lindroos <robert.lindroos at ki.se>
# 
# Original version by 
# Alexander Kozlov <akozlov at kth.se>
# Kai Du <kai.du at ki.se>
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
import glob
import matplotlib.pyplot as plt



mod        = "./mod/"
defparams  = "./params-rob.json"
morphology = "./morphology/"



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
    '''
    Used for setting the maximal conductance of a segment.
    Scales the maximal conductance based on somatic distance and distribution type.
    
    Parameters:
    d3   = distribution type:
         0 linear, 
         1 sigmoidal, 
         2 exponential
         3 step function
    dist = somatic distance of segment
    a4-7 = distribution parameters 
    g8   = maximal conductance of channel
    
    '''
    
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
            
            
            
def calc_rand_Modulation(mod_list, range_list=False, distribution='uniform'):
    
    '''
    uses numpy to draws random modulation factors in range [0,2] from a given distribution:
        uniform          (same probability for all values)
        centered         (highest probability for values at center of interval)
        inverse centered (higest probability for values close to limits)
    for each channel in mod_list.
    The factors can also be mapped to an arbitrary interval. 
    This is done if a range_list is given.
    
    mod_list     = list of channels to be modulated
    range_list   = list of [min, max] values to be used in modulation
    distribution = distribution to draw factors from
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
    
    


def save_obj(obj, name ):
    with open('save/'+ name + '.pkl', 'wb') as f:
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
                            chan2mod=['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can'] ): 
    
    
    
    
    # initiate cell
    cell = MSN(params=par, factors=factors)
    
    # set cascade ---- move to MSN def?
    casc = h.D1_reduced_cascade2_0(0.5, sec=cell.soma) # other cascades also possible...
    
    pointer = casc._ref_Target1p       #totalActivePKA    (if full cascade used)
    
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
    
    pka = h.Vector()
    pka.record(pointer)
    
    # peak n dipp parameters
    da_peak   = 1500   # concentration [nM]
    da_tstart = 500    # stimulation time [ms]
    da_tau    = 500    # time constant [ms]
    
    
    tstop = simDur               # [ms]
    
    
    # all channels to modulate
    mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can' ]
    
    
    not2mod = [] #['kaf']
    
    
    # find channels that should not be modulated
    for chan in mod_list:
        
        if chan not in chan2mod:
            
            not2mod.append(chan)
    
    
    # calc modulation factors--------------------------------------------------------------
    base_mod    = casc.init_Target1p
    
    # for random modulation: modValues = np.arange(0.1, 2.0, 0.1) 
    if randMod == 1:
    
        if amp == 0.32:
        
            mod_fact = calc_rand_Modulation(mod_list, range_list=[[0.60,0.80],    \
                                                                  [0.65,0.85],  \
                                                                  [0.75,0.85],  \
                                                                  [0.85,1.25],  \
                                                                  [1.0,2.0],    \
                                                                  [1.0,2.0],    \
                                                                  [0.0,1.0]],
                                                                  distribution='uniform'  )
            
        else:
            
            mod_fact = RES[run]['factors']
    
    else:
        mod_fact = [ 0.8, 0.8, 0.8, 1.25, 2.0, 2.0, 0.5  ]
           
    
    factors     = []
    for i,mech in enumerate(mod_list):
        
        # normalization of factors to substrate range. Only for dynamical mod?
        factor  = (mod_fact[i] - 1) / (2317.1 - base_mod) #2317.1
        
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
            
    
    # dynamical modulation ---------------------------------------------------------------
    if sim == 'modulation':
        
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
        
        #print('inne ', sim)
    
        for sec in h.allsec():
            
            for seg in sec:
                
                for mech in seg:
                
                    if mech.name() in mod_list: 
                        
                        if sec.name().find('axon') < 0: # if 0: no axon modulated; if 10: all sections modulated
                        
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
    
    
    
    # set ampa and nmda epsp's--what compartments to use?
    elif sim == 'plateau':
        
        print('inne ', sim)
        
        dend_name = 'dend[' + str(int(section)) + ']'
        
        for sec in h.allsec():
            
            if sec.name() == dend_name:
                
                x = 0.5
                
                ampa = h.ampa(x, sec=sec)
                ampa.onset = 100
                ampa.gmax  = 5e-3
                #h.setpointer(pointer, 'pka', seg.kaf )
                
                nmda = h.nmda(x, sec=sec)
                nmda.onset = 100
                nmda.gmax  = 10e-2
                #h.setpointer(pointer, 'pka', seg.kaf )
                
                vmL = h.Vector()
                vmL.record(sec(x)._ref_v)
                
                d2soma = int(h.distance(x, sec=sec))
                
                #print(sec.name(), h.distance(seg.x, sec=sec))
    
            
    
    # record Ca dynamics!
    elif sim == 'ca':
        
        print('inne ', sim)
        
        for i,sec in enumerate(h.allsec()):
            
            if sec.name().find('axon') < 0: # don't record in axon
            
                for j,seg in enumerate(sec):
                    
                    sName = sec.name().split('[')[0]
                    
                    # N, P/Q, R Ca pool
                    cmd = 'ca_%s%s_%s = h.Vector()' % (sName, str(i), str(j))
                    exec(cmd)
                    
                    cmd = 'ca_%s%s_%s.record(seg._ref_cai)' % (sName, str(i), str(j))
                    exec(cmd)   
                    
                    # the L-type Ca
                    cmd = 'cal_%s%s_%s = h.Vector()' % (sName, str(i), str(j))
                    exec(cmd)
                    
                    cmd = 'cal_%s%s_%s.record(seg._ref_cali)' % (sName, str(i), str(j))
                    exec(cmd)   
                    
                    # uncomment here if testing kaf blocking effect on bAP
                    #gbar = seg.kaf.gbar
                    #seg.kaf.gbar = 0.8 * gbar
    
    
              
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
    if sim in ['vm', 'directMod', 'modulation']:
        
        '''s = ''
        for chan in not2mod:
            s = chan + s'''
        
        save_vector(tm, vm, ''.join(['./vm_', sim, '_', str(int(amp*1e3)), '.out']) )
        '''
        spikes      = getSpikedata_x_y(tm,vm) 
        amp         = int(amp*1e3) 
        
        
        if amp == 320:
            
            RES[run]                = {}
            RES[run]['factors']     = mod_fact
            
            if run == 0:
            
                RES['channels'] = mod_list
        
        RES[run][amp]   = {'spikes': spikes}'''
        
        
    elif sim == 'plateau':
    
        save_vector(tm, vm, ''.join(['./vm_', sim, str(d2soma), '_dend', str(int(section)), '.out']) )
        save_vector(tm, vmL, ''.join(['./vmL_', sim, str(d2soma), '_dend', str(int(section)), '.out']) )
    
    
    elif sim == 'ca':
        
        # vm
        save_vector(tm, vm, ''.join(['./vm_', sim, '_', str(int(amp*1e3)), '.out']) )
        
        # Ca
        for i,sec in enumerate(h.allsec()):
        
            if sec.name().find('axon') < 0:
            
                for j,seg in enumerate(sec):
                    
                    sName   = sec.name().split('[')[0]
                    
                    vName   = 'ca_%s%s_%s' %  (sName, str(i), str(j)    )
                    v2Name   = 'cal_%s%s_%s' %  (sName, str(i), str(j)    )
                    fName   = 'Ca/Org/ca_%s_%s.out'  %  (str(int(np.round(h.distance(seg.x)))), vName )
                    
                    cmd     = 'save_vector(tm, np.add(%s, %s), %s)' % (vName, v2Name, 'fName' ) 
                    
                    exec(cmd)
                



# if run from terminal...   ===============================================================
if __name__ == "__main__":
    
    sys.argv = ['a', 'ca']
    print('starting sim')
    
    arguments = None
    
    # define modulation vectors--used to investigate m and h gate differences, i.e. shift gates
    '''
    kaf_m_vhalf = [8, 9, 10, 11, 12]
    kaf_h_vhalf = [73.6, 74.6, 75.6, 76.6, 77.6]
    kaf_m_slope = [26, 26.5, 27, 27.5, 28]
    kaf_h_slope = [-8, -8.5, -9, -9.5, -10]
    
    # 
    naf_m_vhalf = [23, 24, 25, 26, 27]
    naf_h_vhalf = [60, 61, 62, 63, 64]
    naf_m_slope = [8.2, 8.7, 9.2, 9.7, 10.2]
    naf_h_slope = [-5, -5.5, -6, -6.5, -7]
    '''
    
    # tau sim 1
    #taum = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15]
    #tauh = [0.14, 0.19, 0.24, 0.29, 0.34, 0.39, 0.44, 0.49, 0.54]
    
    # tau sim 2 vhalf
    #taum = [52, 54, 56, 58, 60, 62, 64]
    #tauh = [26, 28, 30, 32, 34, 36, 38]
    
    # tau sim 3 slope
    #taum = [4, 6, 8, 10, 12]
    #taun = [-23, -28, -33, -38, -43]
    #tauh = [2.5, 3.5, 4.5, 6.5, 7.5]
    
    factors = [1] #list(itertools.product(taum, tauh, taun))
    
    if sys.argv[1] == 'ca':     # single run
    
        current = 2000
        main( par="./params-rob.json",          \
                    amp=current*1e-3,           \
                    modulation=0,               \
                    simDur=200,                 \
                    stimDur=2,                  \
                    sim='ca'                    )
                    # dur increased from 2 to 25
                    # amp decreased from 2000 to 1000
                                                    
                                                    
                                                    
                                                    
    elif sys.argv[1] in ['vm', 'directMod']:  #for current in currents:
    
        mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can']
        
        randMod = 0
        
        currents = np.arange(310,335,10)
        
        RES = {}
        
        #for n in range(400):
        
        files = glob.glob('vm*')
        
        plot = 1
        if plot:
            for f in files:
                
                plt.plot(*np.loadtxt(f, unpack=True))
                
            plt.show()
        else:
        
            for current in currents:
                main( par="./params-rob.json",      \
                        amp=current*1e-3,           \
                        run=1,                      \
                        modulation=0,               \
                        simDur=1000,                \
                        stimDur=900,                \
                        sim=sys.argv[1],            \
                        randMod=randMod,            \
                        chan2mod=mod_list           )
                
                
                if randMod == 1 and n % 5 == 0:
                    
                    print('inside save loop', n)
                
                    if n == 0:
                        
                        mod_fact = RES[0]['factors']
                    
                        ID = ''
                        
                        for i,mech in enumerate(mod_list):
                            
                            ID = ID + mech + str( int(mod_fact[i]*100) )
                
                    save_obj(RES, ''.join(['./StatUnikaf-noAxon-', ID]) )
                
                
            
                        
                        
                        
                        
                        
                        
    elif sys.argv[1] == 'plateau':  #for sec in sections (dendritic):
        
        sections = np.arange(1,57) #[4, 8, 18, 40, 57]
        for sec in sections:
            main( par="./params-rob.json",      \
                    amp=0,                      \
                    modulation=0,               \
                    simDur=500,                 \
                    sim=sys.argv[1],            \
                    stimDur=800,                \
                    section=sec                 )
                        
                        
                        
                        
    elif sys.argv[1] == 'modulation':     # multiple currents/modulations?
        
        mod_list = ['naf', 'kas', 'kaf', 'kir', 'cal12', 'cal13', 'can']
        
        currents = np.arange(300,355,10)
        
        randMod = 1
        
        for n in range(80):
            
            RES = {}
        
            for i,current in enumerate(currents):
                 
                main( par="./params-rob.json",      \
                        amp=current*1e-3,           \
                        run=i,                      \
                        modulation=1,               \
                        simDur=3000,                \
                        stimDur=3000,               \
                        sim=sys.argv[1],            \
                        randMod=randMod,            \
                        chan2mod=mod_list           )
                                
            '''
                 chan2mod=[item for item in mod_list if not item == chan ] \

            '''
        
            if randMod == 1:
            
                mod_fact = RES['factors']
                
                ID = ''
                
                for i,mech in enumerate(mod_list):
                    
                    ID = ID + mech + str( int(mod_fact[i]*100) )
            
                save_obj(RES, ''.join(['RMD-', ID]) )
                                                    
    
                                                    
                                                    
                                                    
                                                    
    
    
    
          
    
        

