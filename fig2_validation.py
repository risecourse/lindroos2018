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
from joblib import Parallel, delayed
import multiprocessing
import numpy                as np
import matplotlib.pyplot    as plt
import plot_functions       as fun
import MSN_builder          as build
import json

import os

if not os.path.exists('Results/FI'):
    os.makedirs('Results/FI')

if not os.path.exists('Results/Ca'):
    os.makedirs('Results/FI')
    
h.load_file('stdlib.hoc')
h.load_file('import3d.hoc')

usepar=0 # 0: run in serial (slow, but no pickle errors)
         # 1: run in parallel (fast but requires pickling the tasks)
    
def save_vector(t, v, outfile):
    
    with open(outfile, "w") as out:
        for time, y in zip(t, v):
            out.write("%g %g\n" % (time, y))                     
 
                  




def main(   par="./params_dMSN.json",        \
                            sim='vm',       \
                            amp=0.265,      \
                            run=None,       \
                            simDur=1000,    \
                                parfile=None, \
                            stimDur=900     ): 
    
    
    # initiate cell
    cell    =   build.MSN(  params=par,                             \
                            par = parfile, \
                            morphology='latest_WT-P270-20-14ak.swc' )
    
    
    # set cascade--not activated in this script, 
    # but used for setting pointers needed in the channel mechnisms
    casc    =   h.D1_reduced_cascade2_0(0.5, sec=cell.soma) 
    
    
    # set pointer target in cascade
    pointer =   casc._ref_Target1p    
       
    
    # set edge of soma as reference for dendritic distance 
    h.distance(1, sec=h.soma[0])
    
    
    # set current injection
    stim        =   h.IClamp(0.5, sec=cell.soma)
    stim.amp    =   amp  
    stim.delay  =   100
    stim.dur    =   stimDur    
     
    
    # record vectors
    tm  = h.Vector()
    tm.record(h._ref_t)
    vm  = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)
    
    tstop       = simDur 
    # dt = default value; 0.025 ms (25 us)
                  
    
    # set pointers; need since same mechanisms are used for dynamic modulation of channels.
    # Modulation of channels is not used in this script
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
                
                if sec.name().find('soma') >= 0:
                    
                    
                    # N-type Ca (can) is only distributed to the soma section
                    h.setpointer(pointer, 'pka', seg.can )
                    
            
    
    # configure simulation to record from both calcium pools.
    # the concentration is here summed, instead of averaged. 
    # This doesn't matter for the validation fig, since relative concentration is reported.
    # For Fig 5B, where concentration is reported, this is fixed when plotting.
    # -> see the plot_Ca_updated function in plot_functions.
    if sim == 'ca':
        
        print('configure', sim, 'simulation')
        
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
                    #block_fraction = 0.2
                    #gbar           = seg.kaf.gbar
                    #seg.kaf.gbar   = (1 - block_fraction) * gbar
    
    
              
    # solver------------------------------------------------------------------------------            
    cvode = h.CVode()
    
    h.finitialize(cell.v_init)
    
    # run simulation
    while h.t < tstop:
                
        h.fadvance()
        
    
    # save output ------------------------------------------------------------------------
    
    if sim == 'ca':
        
        print('saving', sim, 'simulation')
        
        # vm
        save_vector(tm, vm, ''.join(['Results/Ca/vm_', sim, '_', str(int(amp*1e3)), '.out']) )        
        
        # ca
        for i,sec in enumerate(h.allsec()):
        
            if sec.name().find('axon') < 0:
            
                for j,seg in enumerate(sec):
                    
                    
                    sName       =   sec.name().split('[')[0]
                    vName       =   'ca_%s%s_%s'  %  ( sName, str(i), str(j)  )
                    v2Name      =   'cal_%s%s_%s' %  ( sName, str(i), str(j)  )
                    fName       =   'Results/Ca/ca_%s_%s.out'  %  ( str(int(np.round(h.distance(seg.x)))), vName )
                    
                    cmd     = 'save_vector(tm, np.add(%s, %s), %s)' % (vName, v2Name, 'fName' ) # this is were concentrations are summed (see above)
                    
                    exec(cmd)
        
                    
    elif sim == 'vm':
        
        print('saving', sim, 'simulation', str(int(amp*1e3)))
        
        # vm
        save_vector(tm, vm, ''.join(['Results/FI/vm_', sim, '_', str(int(amp*1e3)), '.out']) )
                


# Start the simulation.
# Function needed for HBP compability  ===================================================
if __name__ == "__main__":
    
    
    print('starting sim')
    par="./params_dMSN.json"
    with open(par) as file:
        parfile = json.load(file)
    
    # dendritic validation: change in [Ca] following a bAP (validated against Day et al., 2008)
    current = 2000
    main( par="./params_dMSN.json",          \
         parfile = parfile, \
                amp=current*1e-3,           \
                simDur=200,                 \
                stimDur=2,                  \
                sim='ca'                    )
                                                
    
    print('starting somatic excitability simulation')                                               
    if usepar==1:
        # somatic excitability (validated against FI curves in Planert et al., 2013)  
        currents    = np.arange(-100,445,40)
        num_cores   = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(main)(   par="./params_dMSN.json",   \
                                                    parfile = parfile, \
                                                    amp=current*1e-3,           \
                                                    run=1,                      \
                                                    simDur=1000,                \
                                                    stimDur=900                 \
                            ) for current in currents)
                            
        currents    = np.arange(320,445,40)
        Parallel(n_jobs=num_cores)(delayed(main)(   par="./params_dMSN.json",   \
                                                    parfile = parfile, \
                                                    amp=current*1e-3,           \
                                                    run=1,                      \
                                                    simDur=1000,                \
                                                    stimDur=900                 \
                            ) for current in currents)
    else:
    # somatic excitability (validated against FI curves in Planert et al., 2013)  
        currents    = np.arange(-100,445,40)
        for current in currents:
            main(   par="./params_dMSN.json",   \
                  parfile = parfile, \
                   amp=current*1e-3,           \
                    run=1,                      \
                    simDur=1000,                \
                    stimDur=900 
                    )
                            
        currents    = np.arange(320,445,40)
        for current in currents:
            main(   par="./params_dMSN.json",   \
                    parfile = parfile, \
                    amp=current*1e-3,           \
                    run=1,                      \
                    simDur=1000,                \
                    stimDur=900                 \
                    )

                        
    
    print('all simulations done! Now plotting')
        
    # PLOTTING
    fun.plot_Ca('Results/Ca/ca*.out')
    fun.plot_vm()
    plt.show()        

                                                    
    
                                                    
                                                    
                                                    
                                                    
    
    
    
          
    
        

