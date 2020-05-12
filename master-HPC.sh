#!/bin/sh
# This is the main script used for design an inlet with large gas bubbles entering a continuous liquid. This script is only used to concentrate the user-input for the model (and is dependent on the flow solver the user wants to use, denoted by the main IF-statement in the script. This script was initially developed for an OpenFOAM-case. The modelling itself is done in two Python-scripts: the first one ("TubeBundle_ReadInlet_<solver>.py") is flow solver-dependent and is used to create a table from the inlet geometry of the case. The second script ("TubeBundle_InletModelling.py") is solver-independent and creates the actual inlet. The bubble shapes which are used, are defined in the "TubeBundle-InletModelling" script. 

# User input 
export CFD_PROGRAMME=OpenFOAM #What CFD-solver will be used? Only OpenFOAM is modelled.
export CFD_VERSION=4.1-intel-2017a #Version of the CFD-solver (module name= CFD_PROGRAMME/CFD_VERSION)
export PYTHON_VERSION=Anaconda2/4.2.0 # What kind of Python should be used?
export DIM=3 # Number of geometrical dimensions of the case (2 or 3)
export CASE_PATH=$VSC_SCRATCH/TB-timing #Location of base case
export startTime=0 #Time step folder from where inlet should be modelled
export endTime=1 #Last flow time to be defined
export timeStepSize=0.001 #Time step size to be used in following calculation
export tunit=0.25 #Unit time scale - parameter for inlet modelling
export inletName=inlet #Name of the inlet boundary to be modelled
export rhog=1.2 #Density of the gas
export rhol=1000 #Density of the liquid
export mg=0.00005 #Amount of the gas to be introduced in domain over a time 'tunit'
export tol_mg=1e-07 #Tolerenace on the amount of gas to be introduced in domain over a time 'tunit' 
export U=1.5 #Velocity of the mixture  to be introduced in domain 
export intersectBoundary=True #Boolean indicating whether prescribed bubble shapes may intersect with the domain boundaries
export intersectBubble=False #Boolean indicating whether prescribed bubble shapes may intersect with previously defined bubbles

# Load Python (inlet modelling performed in Python Anaconda)
module load $PYTHON_VERSION

# Execution of the script depends on CFD-programme to be used
export CFD_MODULE=$CFD_PROGRAMME/$CFD_VERSION
if [ "$CFD_PROGRAMME" = "OpenFOAM" ]
then
	# Read the inlet geometry : This script is designed to be case-dependent
	python TubeBundle_readInlet_$CFD_PROGRAMME.py $DIM $CASE_PATH $CFD_MODULE $startTime $inletName
else
	echo "The CFD-programme you have tried to use is not defined in the main script."
	exit 1
fi

# Model the inlet : This script is designed to be case-independent
python TubeBundle_inletModelling.py $CASE_PATH $startTime $endTime $timeStepSize $tunit $inletName $rhog $rhol $mg $tol_mg $U $intersectBoundary $intersectBubble

# The modelled inlet has to be written in a format compatible with the flow solver that is to be used
if [ "$CFD_PROGRAMME" = "OpenFOAM" ]
then
        # Write the inlet boundary condition with Python-script
	python TubeBundle_writeBC_$CFD_PROGRAMME.py $CASE_PATH $startTime $inletName
else
        echo "The CFD-programme you have tried to use is not defined in the main script."
        exit 1
fi


