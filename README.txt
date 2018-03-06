Model used in 

Lindroos et al. (2018). doi 10.3389/fncir.2018.00003

Errata:
There is a discrepancy between the model description of the soma in the paper (Morphology;
page 6) and how the morphology file is interpreted by NEURON.
In the paper it is stated that two additional somatic points are added to the morphology,
giving a length and diameter of 12.2 and 11.2 respectively. The extra points are however
not interpreted as part of the same section as the root node, but as an additional somatic
protrusion. This have the consequence that the somatic area, and thereby currents, are 
larger than expected. The total somatic area of the model, instead corresponds to the size 
of a sphere with a radius of about 7 um. 

This does not affect the results of the paper since the same result can be obtained in a 
model with sperical soma of arbitrary size. This is done by rescaling the maximal value of 
the somatic conductances by the following scale factor:
    
    SF = A_org/A_new
    
Where A_org is the original area (626.3078906869523 um2, in this case), and
      A_new is the somatic area of the model you wish to mapp to.
      
Since the dendritic conductances are set based on the distance to the edge of the soma, 
the distance of each segment must also be rescaled based on the difference in somatic 
(root) radius:

    dist(seg_x) = dist + (r_new - r_org)
    
 

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



 
