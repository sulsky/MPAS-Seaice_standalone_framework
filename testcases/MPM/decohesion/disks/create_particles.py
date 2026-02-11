from netCDF4 import Dataset
import numpy as np
import math
import sys
sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions
import os
import argparse

#--------------------------------------------------------------------

def create_particles(gridFilename, meshType):

    if (meshType == 'quad'):
       particleInitType = "number"
       particleInitNumber = "4"
       particlePositionInitType = "even"
       particleGeometry = "disks"
       sphereRadius = 1.0
       filenameOut = "particles.nc"
       initial_particle_positions(gridFilename,
                               filenameOut,
                               particleInitType,
                               particleInitNumber,
                               particlePositionInitType,
                               particleGeometry,
                               sphereRadius)
    elif (meshType == 'hex'):
       particleInitType = "number"
       particleInitNumber = "9"
       particlePositionInitType = "onePerEdge"
       particleGeometry = "disks"
       sphereRadius = 1.0
       filenameOut = "particles.nc"
       initial_particle_positions(gridFilename,
                                filenameOut,
                                particleInitType,
                                particleInitNumber,
                                particlePositionInitType,
                                particleGeometry,
                                sphereRadius)

    # read position data from particle file
    filenameIn = filenameOut
    fileIn = Dataset(filenameIn, "r")
    nParticles = len(fileIn.dimensions["nParticles"])
    posnMP = fileIn.variables["posnMP"][:]
    fileIn.close()

    # create velocity at material points
    uvVelMP = np.zeros([nParticles, 2])

    for iParticle in range(0, nParticles):
           xp = posnMP[iParticle, 0]

           if (xp < 0.5):
               uvVelMP[iParticle, 0] = 0.1
               uvVelMP[iParticle, 1] = 0.1
           else:
               uvVelMP[iParticle, 0] = -0.1
               uvVelMP[iParticle, 1] = -0.1

    # append velocity to particle file
    filenameOut = filenameIn
    fileOut = Dataset(filenameOut, "a", format="NETCDF3_CLASSIC")
    var = fileOut.createVariable("uvVelMP","d",dimensions=["nParticles", "TWO"])
    var[:] = uvVelMP[:]
    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='gridFilename', required=True)
    parser.add_argument('-m', dest='meshType', required=True)

    args = parser.parse_args()

    create_particles(args.gridFilename, args.meshType)
