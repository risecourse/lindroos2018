#!/usr/bin/env python

'''
cleans up meassy mod file.

After exporting cascade models to sbml from matlab (simbiology) or copasi the
resulting parameter and state names are long and and hard to read. When further
converting these files to the neuron readable .mod files (using NeuroML) the 
resulting files are close to impossible to comprehend.

This script maps the ID's from the mod file with NAME's in the xml file making the final 
script human readable.

The same clean up is also performed the nrn.py file resulting from the 
neuroML -> mod conversion 


v.1.2 
written by Robert Lindroos Feb-May. 2016
robert . lindroos at ki . se

'''
import sys
import re

# in file 2 manipulate
#fname =  '../D1_Neuron_Model_0.mod'
fname = sys.argv[1]
py = sys.argv[3]

# dict to store parameters in
parameters = {}


flag = '0'      # flag to tell if we are in NEURON or STATE block
counter = 1     # counter used for setting new parameter/variable names

# function for automatic keys
def chosePname(p):       
    
    global counter
    
    p = p + str(counter)      # new name is p<n>. eg p1, s14
    
    counter += 1        
    return p
    
# function that extracts parameter/species name from parameterfile
def getName(s_id):

    xml = open(sys.argv[2], 'r')
    
    for line in xml.readlines(): # for line in file...
        
        tags = line.split()
        
        if not tags:
            
            continue
        
        if tags[0] in ['<compartment', '<species', '<parameter', '<model']:
            
            Dict = {}
            
            for tag in tags:
                
                split = tag.split('=')
                
                if len(split) > 1:
                    
                    Dict[split[0]] = split[1].strip('"')
            
            #print Dict['id'], s_id    
            if Dict['id'] == s_id:
                xml.close()
                return Dict['name']
    xml.close()                 
    return 0
                        

#### read in parameters to be changed #################################

with open(fname) as f:
    for line in f.readlines():  # for line in file...
            
        l = line.split()
        
        if not l:
            continue
        
        if l[0] == 'NEURON':     # if NEURON or STATE block - start saving parameters
            
            flag = 'N'
            
        elif l[0] == 'STATE':
            
            flag = 'S'
            
        elif flag != '0':
             
            if l[0] == '}':     # end of block--reset to non save mode
                flag = '0'
                counter = 1
                continue
                
            if flag == 'N' and l[1][0:2] == 'mw':  
                
                newParameter = chosePname('p')     # use function to set parameter value    
                parameters[newParameter] = l       # add line list to dict
                
            elif flag == 'S' and l[0][0:2] == 'mw':
                
                newParameter = chosePname('s')     # use function to set state value    
                parameters[newParameter] = l       # add line list to dict
                



#### open file, create back up and execute the changes #################################

# open and read file 
f = open(fname,'r')
filedata = f.read() 

newdata = filedata

# create backup
f2 = open(''.join([fname,'.bak']),'w')
f2.write(filedata)
f2.close() 

# open the nrn.py file
pyf = open(py,'r')
pydata = pyf.read()   

flag = 0
for i,key in enumerate(parameters):
    
    # get ID from dict
    if len(parameters[key]) > 1:    # if parameter from NEURON block
        ID = parameters[key][1]
    else:
        ID = parameters[key][0]    # else--STATE block

    # "clean up" ID    
    if ID[-2:] == '_0':
        ID = ID[:-2]  # remove the last two characters (_0) since don't excits in xml file
        
    elif ID[-12:-2] == '_reaction_':
        ID = ID[:-12]  # remove reaction ID from parameter--only works if reactions are less than 10.
        flag = '_reaction'
        
    # if exists--remove colon from org string
    ID = re.sub('\:$', '', ID) 

    # get names from xml file--error later if not found...   
    NAME = getName(ID)      
    print NAME, ID
    # if exists exchange * with _ in names
    NAME = re.sub('\*', '_', NAME)
    
    if not flag == 0:
        ID = ID + flag
        flag = 0
        
    # search file and replace ID with name
    newdata = newdata.replace(ID,NAME)
    pydata = pydata.replace(ID,NAME)

# replace any trailing _0 with _O
newdata = newdata.replace('_0', '_model') 
pydata = pydata.replace('_0', '_model')  

# shorten long variable names (e.g. rate_irreversable_...)
newdata = newdata.replace('rate__revreaction', 'r_r')
newdata = newdata.replace('rate__irrevreaction', 'r_ir')
#newdata = newdata.replace('rate__reaction', 'r_r')

f.close()
pyf.close()

# re-write to file
f = open(fname,'w')
f.write(newdata)
f.close()

pyf = open(py,'w')
pyf.write(pydata)
pyf.close()

print pydata

print 'replacement finished!'
        
        
        
        
 
        
        
    
    
    
    
        
        
