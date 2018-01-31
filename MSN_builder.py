#
'''
The MSN class defining the cell
'''

from neuron import h
from math import exp
import json

# Distributions:
'''
T-type Ca: g = 1.0/( 1 +np.exp{(x-70)/-4.5} )
naf (den): (0.1 + 0.9/(1 + np.exp((x-60.0)/10.0)))

'''

def calculate_distribution(d3, dist, a4, a5,  a6,  a7, g8):
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

            
        

# ======================= the MSN class ==================================================

class MSN:
    def __init__(self, params=None):
        Import = h.Import3d_SWC_read()
        Import.input('latest_WT-P270-20-14ak.swc')
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
        for sec in self.axonlist:
            sec.insert('naf')
            sec.insert('kas')
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

        self.distribute_channels("soma", "gbar_naf", 0, 1, 0, 0, 0, float(par['gbar_naf_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kaf", 0, 1, 0, 0, 0, float(par['gbar_kaf_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kas", 0, 1, 0, 0, 0, float(par['gbar_kas_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kdr", 0, 1, 0, 0, 0, float(par['gbar_kdr_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kir", 0, 1, 0, 0, 0, float(par['gbar_kir_somatic']['Value']))
        self.distribute_channels("soma", "gbar_sk",  0, 1, 0, 0, 0, float(par['gbar_sk_somatic']['Value']))
        self.distribute_channels("soma", "gbar_bk",  0, 1, 0, 0, 0, float(par['gbar_bk_somatic']['Value']))
        
        self.distribute_channels("axon", "gbar_naf", 3, 1, 1.1, 30, 500, float(par['gbar_naf_axonal']['Value']))
        self.distribute_channels("axon", "gbar_kas", 0, 1, 0, 0, 0,      float(par['gbar_kas_axonal']['Value']))
        
        self.distribute_channels("dend", "gbar_naf", 1, 0.1, 0.9, 60.0, 10.0, float(par['gbar_naf_basal']['Value']))
        self.distribute_channels("dend", "gbar_kaf", 1,   1, 0.5,  120,  -30, float(par['gbar_kaf_basal']['Value']))
        self.distribute_channels("dend", "gbar_kas", 2,   1, 9.0,  0.0, -5.0, float(par['gbar_kas_basal']['Value']))
        self.distribute_channels("dend", "gbar_kdr", 0, 1, 0, 0, 0, float(par['gbar_kdr_basal']['Value']))
        self.distribute_channels("dend", "gbar_kir", 0, 1, 0, 0, 0, float(par['gbar_kir_basal']['Value']))
        self.distribute_channels("dend", "gbar_sk",  0, 1, 0, 0, 0, float(par['gbar_sk_basal']['Value']))
        self.distribute_channels("dend", "gbar_bk",  0, 1, 0, 0, 0, float(par['gbar_bk_basal']['Value']))

        self.distribute_channels("soma", "pbar_cal12", 0, 1, 0, 0, 0, 1e-5)
        self.distribute_channels("soma", "pbar_cal13", 0, 1, 0, 0, 0, 1e-6)
        self.distribute_channels("soma", "pbar_car",   0, 1, 0, 0, 0, 1e-4)
        self.distribute_channels("soma", "pbar_can",   0, 1, 0, 0, 0, 3e-5)
        self.distribute_channels("dend", "pbar_cal12", 0, 1, 0, 0, 0, 1e-5)
        self.distribute_channels("dend", "pbar_cal13", 0, 1, 0, 0, 0, 1e-6)
        self.distribute_channels("dend", "pbar_car",   0, 1, 0, 0, 0, 1e-4)
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
        '''find maximal dendritic branch length''' 
        h.distance(sec=self.soma)
        dmax = 0
        for sec in self.allseclist:
	        if axon_excluding and sec.name().find('axon') == 0: continue
                dmax = max(dmax, h.distance(1, sec=sec))
        return dmax

    def distribute_channels(self, as1, as2, d3, a4, a5, a6, a7, g8):
        h.distance(sec=self.soma)
        dmax = self._max_dist()
        
        for sec in self.allseclist:
            
            # if right cellular compartment (axon, soma or dend)
            if sec.name().find(as1) >= 0:
                for seg in sec:
                    dist = h.distance(seg.x, sec=sec)
                    val = calculate_distribution(d3, dist, a4, a5, a6, a7, g8)
                    cmd = 'seg.%s = %g' % (as2, val)
                    exec(cmd)
