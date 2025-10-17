# eIF4F Molecular Dynamics Simulation

# Version
numpy 1.26.4
gsd 2.9.0
hoomd 3.5.0

## Set up parameter
The parameter space is included in "parentparam.yml", the parameter could be either ON/OFF, or string or a list of values.

## Generate folder with given parameters
Run
$ py param.py parentparam.yml

It will iterate all parameters and generate folders with a json file included.

## Run simulation
$./runscript.sh "folder name"
