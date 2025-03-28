from netCDF4 import Dataset
import numpy as np
import argparse

#------------------------------------------------------------------

def set_globalID(creationIndexMP,cellIDCreationMP):

    return (creationIndexMP << 32) + cellIDCreationMP

#-------------------------------------------------------------------------------

def check_particles_moved(particleFile0,
                          particleFile1):

    file1 = Dataset(particleFile0,"r")
    file2 = Dataset(particleFile1,"r")

    nParticles1 = len(file1.dimensions["nParticles"])
    nParticles2 = len(file2.dimensions["nParticles"])

    statusMP1 = file1.variables["statusMP"][0,:]
    statusMP2 = file2.variables["statusMP"][0,:]
    nParticlesStatus1 = np.count_nonzero(statusMP1)
    nParticlesStatus2 = np.count_nonzero(statusMP2)

    creationIndexMP1  = file1.variables["creationIndexMP"][0,:]
    cellIDCreationMP1 = file1.variables["cellIDCreationMP"][0,:]

    creationIndexMP2  = file2.variables["creationIndexMP"][0,:]
    cellIDCreationMP2 = file2.variables["cellIDCreationMP"][0,:]

    posnMP0 = file1.variables["posnMP"][:]
    posnMP1 = file2.variables["posnMP"][:]

    file1.close()
    file2.close()

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

    posnMP0 = posnMP0[0,(statusMP1 == 1)][particlesOrder1]
    posnMP1 = posnMP1[0,(statusMP2 == 1)][particlesOrder2]

    if (posnMP0.shape != posnMP1.shape):
        raise Exception("check_particles_moved: posnMP0.shape != posnMP1.shape")

    if (np.array_equal(posnMP0,posnMP1)):
        raise Exception("check_particles_moved: Particles did not apparently move: "+particleFile0+" to "+particleFile1)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--p0', dest='particleFile0', required=True)
    parser.add_argument('--p1', dest='particleFile1', required=True)

    args = parser.parse_args()

    check_particles_moved(args.particleFile0,
                          args.particleFile1)
