from __future__ import print_function
import argparse
from netCDF4 import Dataset
import numpy as np
import os, sys, os.path, math

#------------------------------------------------------------------

def set_globalID(creationIndexMP,cellIDCreationMP):

    return (creationIndexMP << 32) + cellIDCreationMP

#------------------------------------------------------------------

def add_variable_to_diag_file(file1,variableArray1,variableArray2,variableName,nParticles=None):

    filenameDiag = "vars_differ.nc"

    if (not os.path.isfile(filenameDiag)):
        fileDiag = Dataset(filenameDiag,"w",format="NETCDF3_CLASSIC")
    else:
        fileDiag = Dataset(filenameDiag,"a")

    varIn = file1[variableName]
    for dimension in varIn.dimensions:
        if (dimension not in fileDiag.dimensions):
            if (dimension == "nParticles"):
                fileDiag.createDimension(dimension,nParticles)
            else:
                fileDiag.createDimension(dimension,len(file1.dimensions[dimension]))

    varOut = fileDiag.createVariable(varIn.name, varIn.dtype, varIn.dimensions)
    varOut[:] = variableArray2[:] - variableArray1[:]

    fileDiag.close()

#------------------------------------------------------------------

def compare_variable(variableArray1,variableArray2):
    # return true if same

    return np.array_equal(variableArray1,variableArray2)

#------------------------------------------------------------------

def compare_particle_variable(variableArray1,variableArray2,statusMP1,statusMP2,particlesOrder1,particlesOrder2):
    # return true if same

    if (variableArray1.ndim == 2):
        arr1 = variableArray1[0,(statusMP1 == 1)][particlesOrder1]
        arr2 = variableArray2[0,(statusMP2 == 1)][particlesOrder2]
    elif (variableArray1.ndim == 3):
        arr1 = variableArray1[0,(statusMP1 == 1),:][particlesOrder1]
        arr2 = variableArray2[0,(statusMP2 == 1),:][particlesOrder2]
    elif (variableArray1.ndim == 4):
        arr1 = variableArray1[0,(statusMP1 == 1),:,:][particlesOrder1]
        arr2 = variableArray2[0,(statusMP2 == 1),:,:][particlesOrder2]
    else:
        raise Exception("Wrong number of dimensions: %i" %(variableArray1.ndim))

    if (arr1.shape != arr2.shape):
        raise Exception("arr1.shape != arr2.shape")

    return np.array_equal(arr1,arr2), arr1, arr2

#------------------------------------------------------------------

def compare_files(filename1, filename2, logfile, variableNamesIgnore=[]):

    # init error numbers
    nErrorsNonArray = 0
    nErrorsArray = 0

    # check files exist
    if (not os.path.exists(filename1)):
        logfile.write("File %s does not exists!\n" %(filename1))
        sys.exit()
    if (not os.path.exists(filename2)):
        logfile.write("File %s does not exists!\n" %(filename2))
        sys.exit()

    # open files
    file1 = Dataset(filename1, "r")
    file2 = Dataset(filename2, "r")

    # dimensions comparison
    dimensionNames1 = set(file1.dimensions.keys())
    dimensionNames2 = set(file2.dimensions.keys())

    # check to see if dimensions are in one and not the other
    dimensionNamesIn1Not2 = dimensionNames1.difference(dimensionNames2)
    dimensionNamesIn2Not1 = dimensionNames2.difference(dimensionNames1)

    for dimensionsName in dimensionNamesIn1Not2:
        logfile.write("Dimension found in file 1 and not file 2: %s\n" %(dimensionsName))
        nErrorsNonArray = nErrorsNonArray + 1

    for dimensionsName in dimensionNamesIn2Not1:
        logfile.write("Dimension found in file 2 and not file 1: %s\n" %(dimensionsName))
        nErrorsNonArray = nErrorsNonArray + 1

    # check dimension sizes
    dimensionsNameIntersection = dimensionNames1.intersection(dimensionNames2)

    for dimensionName in dimensionsNameIntersection:

        if (dimensionName != "nParticles"):

            len1 = len(file1.dimensions[dimensionName])
            len2 = len(file2.dimensions[dimensionName])

            if (len1 != len2):

                logfile.write("Dimension sizes differ: %i v %i for %s\n" %(len1, len2, dimensionName))
                nErrorsNonArray = nErrorsNonArray + 1

    # particles
    hasParticles = False
    if ("nParticles" in dimensionNames1 and
        "nParticles" in dimensionNames2):

        nParticles1 = len(file1.dimensions["nParticles"])
        nParticles2 = len(file2.dimensions["nParticles"])

        statusMP1 = file1.variables["statusMP"][0,:]
        statusMP2 = file2.variables["statusMP"][0,:]
        nParticlesStatus1 = np.count_nonzero(statusMP1)
        nParticlesStatus2 = np.count_nonzero(statusMP2)

        if (nParticlesStatus1 == nParticlesStatus2):
            hasParticles = True

            creationIndexMP1  = file1.variables["creationIndexMP"][0,:]
            cellIDCreationMP1 = file1.variables["cellIDCreationMP"][0,:]

            creationIndexMP2  = file2.variables["creationIndexMP"][0,:]
            cellIDCreationMP2 = file2.variables["cellIDCreationMP"][0,:]

            globalID1 = []
            globalID2 = []
            for iParticle in range(0,nParticles1):
                if (statusMP1[iParticle] == 1):
                    globalID = set_globalID(creationIndexMP1[iParticle],cellIDCreationMP1[iParticle])
                    globalID1.append(globalID)
            for iParticle in range(0,nParticles2):
                if (statusMP2[iParticle] == 1):
                    globalID = set_globalID(creationIndexMP2[iParticle],cellIDCreationMP2[iParticle])
                    globalID2.append(globalID)
            globalID1 = np.array(globalID1)
            globalID2 = np.array(globalID2)
            particlesOrder1 = np.argsort(globalID1)
            particlesOrder2 = np.argsort(globalID2)
        else:
            logfile.write("nParticlesStatus differ: %i %i\n" %(nParticlesStatus1,nParticlesStatus2))
            nErrorsNonArray = nErrorsNonArray + 1

    # variables comparison
    variableNames1 = set(file1.variables.keys())
    variableNames2 = set(file2.variables.keys())

    # check to see if variables are in one and not the other
    variableNamesIn1Not2 = variableNames1.difference(variableNames2)
    variableNamesIn2Not1 = variableNames2.difference(variableNames1)

    for variablesName in variableNamesIn1Not2:
        logfile.write("Variable found in file 1 and not file 2: %s\n" %(variablesName))
        nErrorsNonArray = nErrorsNonArray + 1

    for variablesName in variableNamesIn2Not1:
        logfile.write("Variable found in file 2 and not file 1: %s\n" %(variablesName))
        nErrorsNonArray = nErrorsNonArray + 1

    # check variable dimensions
    variablesNameIntersection = variableNames1.intersection(variableNames2)

    for variableName in variablesNameIntersection:

        dimensions1 = set(file1.variables[variableName].dimensions)
        dimensions2 = set(file2.variables[variableName].dimensions)

        dimensionsIn1Not2 = dimensions1.difference(dimensions2)
        dimensionsIn2Not1 = dimensions2.difference(dimensions1)

        for dimensionName in dimensionsIn1Not2:
            if (dimensionName != "nParticles"):
                logfile.write("Variable dimension found in file 1 and not file 2: %s %s\n" %(variableName,dimensionName))
                nErrorsNonArray = nErrorsNonArray + 1

        for dimensionName in dimensionsIn2Not1:
            if (dimensionName != "nParticles"):
                logfile.write("Variable dimension found in file 2 and not file 1: %s %s\n" %(variableName,dimensionName))
                nErrorsNonArray = nErrorsNonArray + 1

    # check variable contents
    for variableName in variablesNameIntersection:

        variable1 = file1.variables[variableName]
        variable2 = file2.variables[variableName]

        variableArray1 = variable1[:]
        variableArray2 = variable2[:]

        shape1 = np.shape(variableArray1)
        shape2 = np.shape(variableArray2)

        rank1 = len(shape1)
        rank2 = len(shape2)

        # check if ranks are the same
        if (rank1 != rank2):

            logfile.write("Array ranks different: %s %i %i\n" %(variableName, rank1, rank2))
            nErrorsNonArray = nErrorsNonArray + 1

        else:

            arrayOK = True

            # check array sizes
            if ("nParticles" not in variable1.dimensions):
                for iDim in range(0,rank1):

                    if (shape1[iDim] != shape2[iDim]):

                        logfile.write("Array %s dim sizes not same for dimension %i: size1 %i, size2 %i\n" %(variableName, iDim+1, shape1[iDim], shape2[iDim]))
                        nErrorsNonArray = nErrorsNonArray + 1
                        arrayOK = False

            if (arrayOK):

                # compare array values
                if (variableName not in variableNamesIgnore):
                    if ("nParticles" in variable1.dimensions and hasParticles):
                        arraysEqual, arr1, arr2 = compare_particle_variable(variableArray1,variableArray2,statusMP1,statusMP2,particlesOrder1,particlesOrder2)
                        if (not arraysEqual):
                            diff = arr2[:] - arr1[:]
                            minVal = np.amin(diff)
                            maxVal = np.amax(diff)
                            L2Norm = np.linalg.norm(diff)
                            L2ErrorNorm = math.sqrt(np.sum(np.power(diff,2))/np.sum(np.power(arr1,2)))

                            logfile.write("Arrays %s differ! min: %g, max: %g, L2: %g L2 rel: %g\n" %(variableName,minVal,maxVal,L2Norm,L2ErrorNorm))
                            add_variable_to_diag_file(file1,arr1,arr2,variableName,nParticlesStatus1)
                            nErrorsArray = nErrorsArray + 1

                    if ("nParticles" not in variable1.dimensions and
                        not compare_variable(variableArray1,variableArray2)):

                        diff = variableArray2[:] - variableArray1[:]
                        minVal = np.amin(diff)
                        maxVal = np.amax(diff)
                        L2Norm = np.linalg.norm(diff)
                        L2ErrorNorm = math.sqrt(np.sum(np.power(diff,2))/np.sum(np.power(variableArray1,2)))

                        logfile.write("Arrays %s differ! min: %g, max: %g, L2: %g L2 rel: %g\n" %(variableName,minVal,maxVal,L2Norm,L2ErrorNorm))
                        add_variable_to_diag_file(file1,variableArray1,variableArray2,variableName)
                        nErrorsArray = nErrorsArray + 1


    # close files
    file1.close()
    file2.close()

    # return error numbers
    return nErrorsArray, nErrorsNonArray

#------------------------------------------------------------------

# compare_mpas_files.py -f1 testFile1.nc -f2 testFile2.nc
# compare_mpas_files.py -f1 testFile1.nc -f2 testFile2.nc -i ignoreVarnamesTest.txt

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='BFB comparison for MPAS files.')

    parser.add_argument('-f1', help='First filename to compare.',             dest="filename1", required=True)
    parser.add_argument('-f2', help='Second filename to compare.',            dest="filename2", required=True)
    parser.add_argument('-i',  help='Text file of variable names to ignore.', dest="ignoreFile", default=None)

    args = parser.parse_args()

    variableNamesIgnore = []
    if (args.ignoreFile != None):
        fileIgnoreList = open(args.ignoreFile,"r")
        variableNamesIgnore = fileIgnoreList.readlines()
        fileIgnoreList.close()
    variableNamesIgnore = [word.strip() for word in variableNamesIgnore]

    logfile = open("log_test.txt")

    nErrorsArray, nErrorsNonArray = compare_files(args.filename1, args.filename2, logfile, variableNamesIgnore)
    print("Number of array errors:     ", nErrorsArray)
    print("Number of non-array errors: ", nErrorsNonArray)

    logfile.close()
