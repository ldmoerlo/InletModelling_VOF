# This script is developed to read in the output from the flow solver-independent "inletModelling" script., search for a specific inlet boundary and rewrite the cell centres and faces on that boundary to an array to be used for the inlet modelling in a later stage (not in this script)
# This output contains two matrices, one for flow speed and one for Volume-Of-Fluid of water, where the rows are the cell IDs and the columns are the timesteps for which U/VOFw are defined.
# These matrices are converted into text files, to be written in the OpenFOAM-directory such that the OpenFOAM-utility "timeVaryingMappedFixedValue" can be used.
# This script is called automatically by the masterscript 'TubeBundle_master.sh', so the user input is channeled to this python script from the bash-script directly.

# import of utilities
import numpy as np
import sys
import os # to be able to run Linux terminal commands

# Read input from bash-script
if len(sys.argv) != 4:
    sys.exit("3 arguments should be used: case path - folder of first timestep - inlet boundary name.  ")
casePath=str(sys.argv[1])
firstTimeStep=str(sys.argv[2])
inletName=str(sys.argv[3])

#Separate function definitions
#writeHeader: to write OpenFOAM-header in file at location 'fileLoc' - class and object of parameter should be given to function
def writeHeader(fileLoc,className,objectName):
    f=open(fileLoc,'w')
    f.write(r'/*--------------------------------*- C++ -*----------------------------------*\\'+"\n")
    f.write(r'| =========                 |                                                 |'+"\n")
    f.write(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |'+"\n")
    f.write(r'|  \\    /   O peration     | Version:  4.x                                   |'+"\n")
    f.write(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |'+"\n")
    f.write(r'|    \\/     M anipulation  |                                                 |'+"\n")
    f.write(r'\*---------------------------------------------------------------------------*/'+"\n")
    f.write(r'FoamFile'+"\n")
    f.write(r'{'+"\n")
    f.write('\t version \t\t 4.1;'+"\n")
    f.write('\t format \t\t ascii;'+"\n")
    f.write('\t class \t\t ' + className + ';'+"\n")
    f.write('\t object \t\t ' + objectName + ';'+"\n");
    f.write('}'+"\n")
    f.write(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'+"\n")
    f.write("\n")  
    f.close()
    
#writeFooter: to write OpenFOAM-footer in file at location 'fileLoc'
def writeFooter(fileLoc):
    f=open(fileLoc,'a+')
    f.write("\n")
    f.write(r'// ************************************************************************* //'+"\n")
    f.close()

# Inlet model has been previously stored in "$CASE_PATH/inletDefinition-U.npy" and "$CASE_PATH/inletDefinition-VOFw.npy"
# Time instants at which these variables were defined are stored in "$CASE_PATH/inletDefinition-time.npy"
# Coordinates of the inlet geometry points are stored in "$CASE_PATH/inletPython.npy"
print("Reading defined inlet from Python-files")
UVal=np.load(casePath+'/inletDefinition-U.npy') # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of velocity 
VOFwVal=np.load(casePath+'/inletDefinition-VOFw.npy') # Matrix containing 'coordList' rows (#cell centers) and 'timeVal' columns (# time steps defined) - value of VOFw
timeVal=np.load(casePath+'/inletDefinition-time.npy') # List containing the time instants where U and alpha.water are defined
coordList=np.load(casePath+'/inletPython.npy')
print("Finished reading Python-files. \n")
 
 
# # Test arrays 
# coordList=np.array([[0,-1,-1,0,0.1],[1,-1,1,0,0.1],[2,1,1,0,0.1],[3,1,-1,0,0.1]])
# UVal=np.zeros([4,2,3])
# for j in np.arange(len(UVal[:,0,0])):
#     UVal[j,0,0]=0
#     UVal[j,0,1]=0
#     UVal[j,0,2]=1
#     UVal[j,1,0]=0
#     UVal[j,1,1]=0
#     UVal[j,1,2]=2        
# VOFwVal=np.zeros([4,2,1])
# for j in np.arange(len(VOFwVal[:,0,0])):
#     VOFwVal[j,0,0]=1
#     VOFwVal[j,1,0]=0
# timeVal=np.array([0,1])


# This code writes the inlet values assuming the boundary condition 'timeVaryingMappedFixedValue' is used. It is first checked whether this BC is indeed used. Otherwise, the script returns an error.
# Also, it is checked whether an averaging operation ('setAverage    true;' in OF) is defined - if so, another error is returned.
print("Checking inlet definition in folder "+ firstTimeStep +".")
# Check boundary condition for 'U'
os.system("cd " + casePath + "; grep -nr " + inletName + " "+ firstTimeStep + "/U | cut -d : -f 1 > lineNr_U" ) 
lineNameNr=int(open(casePath+"/lineNr_U",'r').readline())
lineTypeU=lineNameNr+1 # inlet BC type is defined on this line - considering Python starts at index zero
readTypeU=(open(casePath+"/"+firstTimeStep+"/U").readlines())[lineTypeU]
boundaryConditionU=(readTypeU.split())[-1][0:-1]
os.system("rm " + casePath +"/lineNr_U") 
if boundaryConditionU != "timeVaryingMappedFixedValue":
    sys.exit("The condition for 'U' at the boundary '" + inletName +"' is not set to 'timeVaryingMappedFixedValue'. \n")
os.system("cd " + casePath + "; grep -nr " + 'setAverage' + " "+ firstTimeStep + "/U | cut -d : -f 1 > lineNr_setAvg" ) 
lineNameNr=int(open(casePath+"/lineNr_setAvg",'r').readline())-1 #Line where 'setAverage' is defined
readSetAvg=(open(casePath+"/"+firstTimeStep+"/U").readlines())[lineNameNr]
setAvg=(readSetAvg.split())[-1][0:-1]
os.system("rm " + casePath +"/lineNr_setAvg")    
if setAvg != "false":
    sys.exit("Error! The boundary condition at '" + inletName + "' defines an averaging operation for 'U'. This is not compatible with the transient inlet modelling defined in the Python script. \n")
# Check boundary condition for 'alpha.water'
os.system("cd " + casePath + "; grep -nr " + inletName + " "+ firstTimeStep + "/alpha.water | cut -d : -f 1 > lineNr_VOFw" ) 
lineNameNr=int(open(casePath+"/lineNr_VOFw",'r').readline())
lineTypeVOFw=lineNameNr+1 # inlet BC type is defined on this line - considering Python starts at index zero
readTypeVOFw=(open(casePath+"/"+firstTimeStep+"/alpha.water").readlines())[lineTypeVOFw]
boundaryConditionVOFw=(readTypeVOFw.split())[-1][0:-1]
os.system("rm " + casePath +"/lineNr_VOFw")    
if boundaryConditionVOFw != "timeVaryingMappedFixedValue":
    sys.exit("The condition for 'alpha.water' at the boundary '" + inletName + "' is not set to 'timeVaryingMappedFixedValue'. \n")
os.system("cd " + casePath + "; grep -nr " + 'setAverage' + " "+ firstTimeStep + "/alpha.water | cut -d : -f 1 > lineNr_setAvg" ) 
lineNameNr=int(open(casePath+"/lineNr_setAvg",'r').readline())-1 #Line where 'setAverage' is defined
readSetAvg=(open(casePath+"/"+firstTimeStep+"/alpha.water").readlines())[lineNameNr]
setAvg=(readSetAvg.split())[-1][0:-1]
os.system("rm " + casePath +"/lineNr_setAvg")    
if setAvg != "false":
    sys.exit("Error! The boundary condition at '" + inletName + "' defines an averaging operation for 'alpha.water'. This is not compatible with the transient inlet modelling defined in the Python script. \n")
print("Inlet definition in folder "+ firstTimeStep + " is OK. \n")


# Prepare OpenFOAM-directory
print("Creating folder 'boundaryData'.")
errVal=0 # Integer denoting whether os.system has error (>0: at least one error)
try:
    errVal+=os.system("cd " + casePath + r"/constant; mkdir boundaryData; mkdir boundaryData/"+inletName)
    errVal+=os.system("touch " + casePath + r"/constant/boundaryData/"+inletName+r"/points")
    for i in np.arange(len(timeVal)):
        errVal+=os.system("cd " + casePath + r"/constant/boundaryData/" + inletName +"; mkdir "+ str(timeVal[i]) +";")
        errVal+=os.system("touch " + casePath + r"/constant/boundaryData/"+inletName+"/"+str(timeVal[i])+"/U;")
        errVal+=os.system("touch " + casePath + r"/constant/boundaryData/"+inletName+"/"+str(timeVal[i])+"/alpha.water;")
    if errVal>0:
        sys.exit("The Python code encountered a problem preparing the OpenFOAM-directory. Make sure that the 'constant' directory exists but does not contain a folder 'boundaryData'.")
except: 
    sys.exit("The Python code encountered a problem preparing the OpenFOAM-directory. Make sure that the 'constant' directory exists but does not contain a folder 'boundaryData'.")
print("Folder 'boundaryData' was successfully created. \n")


# Write the values of flow speed and Volume-Of-Fluid of water per timestep in the newly created boundaryData folder
# Do not forget to also create boundaryData/inlet/points, containing the coordinates of the inlet points (location where the speed of VOF has to be applied)
print("Writing inlet definition to folder 'boundaryData'.")
nPoints=len(coordList[:,0])
for i in np.arange(len(timeVal)):
    # Set location of files
    fileLoc_points=casePath+r"/constant/boundaryData/"+inletName+"/points"
    fileLoc_U=casePath + r"/constant/boundaryData/"+inletName+"/"+str(timeVal[i])+"/U"
    fileLoc_VOFw=casePath + r"/constant/boundaryData/"+inletName+"/"+str(timeVal[i])+"/alpha.water"
    # Write OpenFOAM-header
    writeHeader(fileLoc_points,"vectorField","points")
    writeHeader(fileLoc_U,"vectorAverageField","values")
    writeHeader(fileLoc_VOFw,"scalarAverageField","values")
    # Write 'points'-file
    f=open(fileLoc_points,'a+')
    f.write(str(nPoints)+'\n')
    f.write('('+'\n')
    for j in np.arange(nPoints):
        f.write('('+str(coordList[j,1])+' '+str(coordList[j,2])+' '+str(coordList[j,3])+') \n')
    f.write(')'+'\n')
    f.close() 
    # Write 'U'-file
    f=open(fileLoc_U,'a+')
    f.write('//Average'+'\n')
    f.write('(0 0 0)'+'\n'+'\n') # necessary for vectorAverageField, but not used if 'setAverage' is set to false (checked earlier in this script)
    f.write('//Data points'+'\n')
    f.write(str(nPoints)+'\n'+'('+'\n')
    for j in np.arange(nPoints):
        f.write('('+str(UVal[j,i,0])+' '+str(UVal[j,i,1])+' '+str(UVal[j,i,2])+') \n')
    f.write(')'+'\n')  
    f.close() 
    # Write 'alpha.water'-file
    f=open(fileLoc_VOFw,'a+')
    f.write('//Average'+'\n')
    f.write('0'+'\n'+'\n') # necessary for vectorAverageField, but not used if 'setAverage' is set to false (checked earlier in this script)
    f.write('//Data points'+'\n')
    f.write(str(nPoints)+'\n'+'('+'\n')
    for j in np.arange(nPoints):
        f.write(str(VOFwVal[j,i,0])+'\n')
    f.write(')'+'\n')     
    f.close() 
    # Write OpenFOAM-footers
    writeFooter(fileLoc_points)
    writeFooter(fileLoc_U)
    writeFooter(fileLoc_VOFw)

print("Finished writing boundary condition to folder 'boundaryData'. \n")


print("Script 'writeBC_OpenFOAM' completed. \n")
