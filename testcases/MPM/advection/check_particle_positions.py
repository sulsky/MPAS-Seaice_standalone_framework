from netCDF4 import Dataset
import sys
import math

#-------------------------------------------------------------------------------

def check_particle_positions():

    filein = Dataset("./output/output.2000.nc","r")

    Time = len(filein.dimensions["Time"])
    nParticles = len(filein.dimensions["nParticles"])

    statusMP1 = filein.variables["statusMP"][0,:]
    cellIDCreationMP1 = filein.variables["cellIDCreationMP"][0,:]
    creationIndexMP1 = filein.variables["creationIndexMP"][0,:]
    posnMP1 = filein.variables["posnMP"][0,:,:]

    errors = []

    for iTime in range(0, Time):

        statusMP2 = filein.variables["statusMP"][iTime,:]
        cellIDCreationMP2 = filein.variables["cellIDCreationMP"][iTime,:]
        creationIndexMP2 = filein.variables["creationIndexMP"][iTime,:]
        posnMP2 = filein.variables["posnMP"][iTime,:,:]


        particleIDToIndex = {}
        for iParticle2 in range(0, nParticles):
            if (statusMP2[iParticle2] == 1):
                particleIDToIndex[(cellIDCreationMP2[iParticle2],creationIndexMP2[iParticle2])] = iParticle2


        maxDistance = -sys.float_info.max
        for iParticle1 in range(0,nParticles):
            if (statusMP1[iParticle1] == 1):
                iParticle2 = particleIDToIndex[(cellIDCreationMP1[iParticle1],creationIndexMP1[iParticle1])]
                distance = math.sqrt(math.pow(posnMP2[iParticle2,0]-posnMP1[iParticle2,0],2) +
                                     math.pow(posnMP2[iParticle2,1]-posnMP1[iParticle2,1],2) +
                                     math.pow(posnMP2[iParticle2,2]-posnMP1[iParticle2,2],2))
                maxDistance = max(maxDistance,distance)

        print("max distance: ", iTime, maxDistance)

        errors.append(maxDistance)


    filein.close()

    print("Final max distance: ", errors[-1])

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    check_particle_positions()
