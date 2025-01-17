from netCDF4 import Dataset
import numpy as np
import math
import sys
from create_ic import latlon_from_xyz, grid_rotation_forward, velocities_analytical
sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions
import os

#--------------------------------------------------------------------

def create_particles():

    rotateCartesianGrid = True
    r = 1.0

    reses = ["2562","10242","40962","163842"]

    icTypes = ["uniform"]

    tests = ["1-x", "y-z", "latlon", "nonlin"]
    #tests = ["1-x"]

    rotateCartesianGrid = True
    r = 1.0

    for test in tests:
        if (test == "1-x"):
            c = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
        elif (test == "y-z"):
            c = [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0]
        elif (test == "latlon"):
           c = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        elif (test == "nonlin"):
           c = [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0]

        for icType in icTypes:

            print("Particle icType: ", icType, "testcase: ", test)

            for res in reses:

                print("  Res: ", res)

                filenameOut = "particles_%s_%s.nc" %(test,res)

                if (not os.path.isfile(filenameOut)):

                    filenameMesh = "grid.%s.nc" %(res)

                    initial_particle_positions(filenameMesh,
                                           filenameOut,
                                           "number",
                                           9,
                                           "onePerEdge",
                                           icType,
                                           1.0)

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
                        yp = posnMP[iParticle, 1]
                        zp = posnMP[iParticle, 2]

                        R = math.sqrt(xp*xp + yp*yp + zp*zp)

                        xp = xp*r/R
                        yp = yp*r/R
                        zp = zp*r/R

                        xp, yp, zp = grid_rotation_forward(xp, yp, zp, rotateCartesianGrid)
                        lat, lon = latlon_from_xyz(xp, yp, zp, r)
                        u, v = velocities_analytical(xp, yp, zp, lat, lon, c)
                        uvVelMP[iParticle, 0] = u
                        uvVelMP[iParticle, 1] = v

                    # append velocity to particle file
                    filenameOut = filenameIn
                    fileOut = Dataset(filenameOut, "a", format="NETCDF3_CLASSIC")
                    var = fileOut.createVariable("uvVelMP","d",dimensions=["nParticles", "TWO"])
                    var[:] = uvVelMP[:]
                    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles()
