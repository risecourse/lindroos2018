Model used in 

Lindroos et al. (2018). doi 10.3389/fncir.2018.00003etan

------------------------------------------------------------------------------------------
ERRATA:

Two corrections were made in the model code to match the description in the article. These 
modifications had minor effects on the results and did not change the conclusions in the 
published paper. The corrected model code closely reproduces the data and figures of the 
article.


1.) The somatic area was larger in the original model than that of a cylinder with length 
    12.2 and diameter 11.2 um (as stated in the paper), but still within reported range 
    (Reinius et al., 2015; Steiner and Tseng, 2017).NEURON interpreted the three point 
    soma as two sections, rather than one. CORRECTIONS MADE: The three point morphology 
    was replaced with a single point soma with radius of 7.06  um, giving the same area as 
    in the original model used in the publication. One additional axonal point was added 
    to the morphology to retain the length of the axon initial segment.

    Additionally, the change in somatic radius also has the consequence that the distance 
    dependence of the dendritic channels are off by 1.46 micrometer compared to the paper 
    version. CORRECTION: The distance dependence of the channels have been recalculated  
    using the following equation:
    
        dist(seg_x) = dist(seg_x) - r_new + r_org = dist(seg_x) - 7.06 + 5.6

2.) A mistake was discovered in the function setting the background synaptic currents that 
    redirected the stimulus intended for the gaba synapse to the glutamate synapse. 
    CORRECTIONS: The pure excitatory drive was replaced with excitation and inhibition in 
    accordance with the paper. Amplitude and activation frequency was rescaled to retain a 
    similar in vivo like state. The rescaling was further done in a way that kept the 
    balance between glutamatergic and gabaergic input stated in the paper.
    
    In general a modification of the GABA/Glut ratio will affect the overall excitability 
    levels, but even with significant changes in this background noise Glut/GABA ratio, 
    the qualitative conclusions of the paper hold.
    
    
Reinius et al., (2015). 
Conditional targeting of medium spiny neurons in the striatal matrix.
Steiner and Tseng., (2017).
Handbook of basal ganglia structure and function. Second edition.

------------------------------------------------------------------------------------------


-Dependencies
-How to run
-File structure
-How to convert the cascade
-Additional comments



DEPENDENCIES
------------------------------------------------------------------------------------------
The model is built in NEURON+Python (v7.4, Hines and Carnevale, 1997)
with the following dependencies:

matplotlib
numpy
__future__
json

+
joblib
multiprocessing
os.path

However, with minor adjustments the model can be run without most of these additional 
libraries.

The plot functions have a lot more dependencies than the simulation files. If you don't 
want to install these libraries you can comment out the plotting command at the end of the
example scripts and plot the data using a software of your own choice. If so, you need to 
move the functions

    getSpikedata_x_y and save_obj

from plot_functions.py to the example script
    
    fig3-4_static_modulation.py





HOW TO RUN THE MODEL
------------------------------------------------------------------------------------------
The following instructions should be valid for Linux and Mac systems.

1. download and decompress the model files
2. compile the mod files by the terminal command "nrnivmodl" (inside the model directory):
    
        -open a terminal
        -cd path/to/directory/
        -nrnivmodl
        
    If your're in the right directory and the your NEURON installation is correct this 
    should give you a lot of text ending with a success message:
        
        Successfully created x86_64/special
        
3. In the same terminal window give the following command:
    
        python <simulation file>
    
    where <simulation file> is one of the three example files described below, e.g. 
        
        fig2_validation.py
        
    This should start the simulation and produce some figure(s) if all dependencies are 
    satisfied. If not install the missing libraries and try again (or restructure the code 
    to remove the dependencies).
    
     



FILE STRUCTURE
------------------------------------------------------------------------------------------
There are three example files that all reproduces some part of the figures in the 
manuscript:

-fig2_validation.py 
Reproduces the model validation in figure 2. It also illustraits how to run the model in
parallel using multiple cores on a local mashine (using the libraries Parallel, delayed 
from joblib and multiprocessing).

-fig3-4_static_modulation.py 
Illustrates how the static DA modulation in figure 3 and 4 were implemented. This 
simulation is implemented to run sequentially and does not depend on the packages above.

-fig6_dynamic_modulation.py
Illustrates how the dynamic modulation was implemented and run (sequentially implemented).

All three depends on the rest of the files to run.

-Ion channels and other mechanisms are stored in the mod files.

-MSN_builder.py 
holds the main class. This is where all active and passive properties are defined and 
channel distributions are set.

-latest_WT-P270-20-14ak.swc 
holds the morphology of the model.

-params-rob.json
holds the parameters used for channel distribution.

-substrates.json
holds initial and maximal values for different substrates in the DA induced intracellular
cascade following a transient DA signal of 500 nM amplitude and 500 ms time constant.
These values are used to linearly map substrate concentration to channel conductance.

plot_functions.py
holds all the functions that were used to create the figures in the manuscript. Only a few
of these are used by the example script


There are also some additional folders:

-Additional_files
holds, for example, the files originally used to run simulations on the Beskow super 
computer. These are a little less well structured than the example files...

-Exp_data
holds experimental data, extracted from published articles (Day et al., 2008 and Planert
et al., 2013), and used for validation.

-Results
these folders are empty. This is where simulation results are stored during simulation.




HOW TO CONVERT THE CASCADE
------------------------------------------------------------------------------------------
The intracellular cascade (Additional_files -> MODEL_speedy_reduced2.xml) can be converted
to a mod file by either one of two strategies, both depending on NeuroML.

1. use the NeuroML tool directly
2. use the python library of NeuroML (pyneuroml)

1.
download/install jNeuroML (https://github.com/NeuroML/jNeuroML).
After successfull installation you should be able to run the following commands (from a
terminal window):
    
    jnml -sbml-import MODEL_speedy_reduced2.xml 1 1
    jnml MODEL_speedy_reduced2_LEMS.xml -neuron
    python sbml_neuroML_mod_parameterCleaner2.py D1_LTP_time_window_0.mod MODEL_speedy_reduced2.xml MODEL_speedy_reduced2_LEMS_nrn.py
    
If neuroml is not globally installed use ./jnml instead of jnml (in local folder)
    
2.
install pyneuroml, lxml and jupyter and open a the notebook TRANSFORM_py_.ipynb 
(Additional_files) from a terminal window:
    
    cd Additional_files/
    jupyter notebook TRANSFORM_py_.ipynb   

and follow the transformation illustrated in the notebook (you might have to move some 
files for the simulation at the end to work).



 
