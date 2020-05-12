# This script is developed to read in an OpenFOAM-case, search for a specific inlet boundary and rewrite the cell centres and faces on that boundary to an array to be used for the inlet modelling in a later stage (not in this script)
# This script is called automatically by the masterscript 'TubeBundle_master.sh', so the user input is channeled to this python script from the bash-script directly.

# import of utilities
import math
import numpy as np
import sys
import os # to be able to run Linux terminal commands
import linecache # read specific line (with known index) from file 
 
 
# Read input from bash-script
if len(sys.argv) != 6:
    sys.exit("5 arguments should be used: dimensions - case path - module of the CFD solver - folder of first timestep - inlet boundary name.  ")
dimensions=str(sys.argv[1])
casePath=str(sys.argv[2])
moduleCFD=str(sys.argv[3])
startTime=str(sys.argv[4])
inletName=str(sys.argv[5])


# Read inlet boundary from OpenFOAM
print("Starting 'writeCellCentres' module of OpenFOAM") 
os.system("cd " + casePath +  "; ml "+ moduleCFD + "; source $FOAM_BASH; writeCellCentres -time " + startTime + ";")
print("Finished 'writeCellCentres'. \n")
print("Loading inlet cell coordinates and face areas into Python.")
for i in np.arange(4):
    if i == 0:
        sourceFile=casePath + "/" + startTime + "/ccx"
    elif i == 1:
        sourceFile=casePath + "/" + startTime + "/ccy"
    elif i == 2:
        sourceFile=casePath + "/" + startTime + "/ccz" 
    elif i == 3:
        sourceFile=casePath + "/" + startTime + "/V" #In the V-file, writeCellCentres writes cell volumes for internal field and patch face areas for boundary fields like the inlet.       
    try:
        # find line where inlet-list starts; you should redo this for all coordinates as 'writeCellCentres' -when possible- replaces lists of uniform values by a single uniform value statement.
        # If the latter is the case, the number of lines you should have in a complete list is obtained from the number of list elements read in another coordinate-file (it's not possible to have >1 'uniform value' coordinate files)
        # The incomplete uniform value-files are subsequently (outside of the try-loop) converted to a list of the same size as the non-uniform list coordinates.
        os.system("cd " + casePath + "; grep -nr " + inletName + " " + sourceFile + " | cut -d : -f 1 > lineNr" ) 
        lineNameNr=int(open(casePath+"/lineNr",'r').readline())
        lineStartNr=lineNameNr+6 # In case of non-uniform list, this is where the list of values in the sourceFile starts
        rowsNrIndex=lineNameNr+4 # On this line, the number of cell centers on the inlet is stated
        os.system("cd " + casePath + "; awk NR==" + str(rowsNrIndex) + " " + sourceFile + " > rowsNr")
        rowsNr=int(open(casePath+"/rowsNr",'r').readline())
        os.system("cd " + casePath + "; rm lineNr rowsNr")
        tempCoordFile=np.ones([rowsNr,1])*float("inf")
        for j in np.arange(rowsNr):
            tempCoordFile[j,0]=float(linecache.getline(sourceFile,lineStartNr+j))
    except ValueError: # If ValueError is triggered, it means that the source-file has a uniform coordinate in the axis you are currently looking
        if not 'rowsNr' in locals(): #If the first coordinate-file has a uniform value, the variable 'rowsNr' does not exist, so you should check whether this variable exists
            checkFile=casePath + "/" + startTime + "/ccy" # if 'rowsNr' does not exist, read the second sourceFile to know the number of rows
            os.system("cd " + casePath + "; grep -nr " + inletName + " " + checkFile + " | cut -d : -f 1 > lineNr" ) 
            lineNameNr_CF=int(open(casePath+"/lineNr",'r').readline())
            rowsNrIndex_CF=lineNameNr_CF+4 # On this line, the number of cell centers on the inlet is stated
            os.system("cd " + casePath + "; awk NR==" + str(rowsNrIndex_CF) + " " + checkFile + " > rowsNr")
            rowsNr=int(open(casePath+"/rowsNr",'r').readline())
            os.system("cd " + casePath + "; rm lineNr rowsNr")
        indexUV=lineNameNr+3
        os.system("cd " + casePath + "; awk NR==" + str(indexUV) + " " + sourceFile + " > unifValue")
        unifValue=float(open(casePath+"/unifValue",'r').readline().split()[-1][0:-1]) #First '-1'm akes sure the value is read, but this still contains a semi-colon, so this should be removed with second index '[0:-1]'.
        os.system("cd " + casePath + "; rm unifValue")
        tempCoordFile=np.ones([rowsNr,1])*float("inf")
        for j in np.arange(rowsNr):
            tempCoordFile[j,0]=unifValue
    if i == 0:
        coordList=np.ones([rowsNr,5])*float("inf") # ID - X - Y - Z - Area
        coordList[:,0]=np.arange(rowsNr)
    coordList[:,(i+1)]=tempCoordFile[:,0]


#Check that all values are inserted correctly
for i in np.arange(rowsNr):
    for j in np.arange(5):
        if coordList[i,j] > 1e15:
            sys.exit("Not all values are correctly read into the Python-script 'TubeBundle_ReadInlet'")
#         echo $((lineNr+rowsNr+5))
#         head -n $((lineNr+rowsNr+5)) $sourceFile > temp
#         tail -n $rowsNr temp > temp2
# done
print("Completed loading of cell coordinates and face areas into Python. \n")


# Determine the normal pointing of the inlet, pointing INTO the domain
# For this purpose, find three linearly independent points in coordList
print("Calculating normal to the inlet. ")
if dimensions == "3":
    tol_product=0.01
    point1=coordList[0,1:4]
    point2=coordList[1,1:4]
    i=2
    point3=coordList[i,1:4]
    while np.linalg.norm(np.cross(point2-point1,point3-point1))/(np.linalg.norm(point2-point1)*np.linalg.norm(point3-point1)) < tol_product:
        i=i+1
        point3=coordList[i,1:4]
    normalInlet=np.cross(point2-point1,point3-point1)*1.0/np.linalg.norm(np.cross(point2-point1,point3-point1))
    # You still need one more point from the domain to determine the correct orientation of the inlet normal
    os.system("cd " + casePath + "; grep -nr '(' constant/polyMesh/points | head -n 1 | cut -d : -f 1 > lineNr" ) 
    lineNameNr=int(open(casePath+"/lineNr",'r').readline())
    lineNr=lineNameNr+1 #  First point that is defined
    os.system("cd " + casePath + "; rm lineNr ")
    i=0
    f=open(casePath+"/constant/polyMesh/points",'r')
    while i < lineNr:
        i=i+1
        f.readline()
    pointDomain=np.double(np.array(f.readline()[1:-2].split(" ")))
    while np.linalg.norm(np.dot(pointDomain-point1,normalInlet))/(np.linalg.norm(pointDomain-point1)) < tol_product: #find pointDomain which is not in the inlet plane - normalInlet already has norm equal to 1
        pointDomain=np.double(np.array(f.readline()[1:-2].split(" ")))
    if np.dot(pointDomain-point1,normalInlet) < 0:
        normalInlet=(-1.0)*normalInlet
    f.close()
elif dimensions == "2":
    normalInlet=np.zeros([3])
    print("The normal to the inlet cannot be calculated directly. Please give the x-, y- and z-coordinates of the normal vector in the following prompts.")
    for i in np.arange(3):
        coordOK=False
        while not(coordOK):
            try:
                temp_coord=float(raw_input("Please provide the "+ str(i+1)+"th coordinate of the normal vector: "))
                coordOK=True
            except ValueError:
                print("Please provide a float value. Try again.")
        normalInlet[i]=temp_coord
    normalInlet=(1/np.linalg.norm(normalInlet))*normalInlet
else :
    sys.exit("Number of dimensions should be either 2 or 3.")
    
print("Completed calculating normal to the inlet, pointing into the domain: "+ str(normalInlet)+".  \n")


# Save inlet and normal in Python Numpy-array format
np.save(casePath+"/inletPython.npy",coordList)
np.save(casePath+"/normalInletPython.npy",normalInlet)


print("Script 'readInlet_OpenFOAM' completed. \n")
