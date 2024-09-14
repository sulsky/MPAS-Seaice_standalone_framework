import sys
sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions

from create_ic import latlon_from_xyz, velocities_strains_analytical,  grid_rotation_forward
from netCDF4 import Dataset
import numpy as np
import math

import os

#--------------------------------------------------------------------

def create_particles():

    mu = 3
    lu = 5

    mv = 2
    lv = 4

    rotateCartesianGrid = True
    r = 1.0

    reses = ["2562","10242","40962","163842"]

    icTypes = ["uniform"]

    for icType in icTypes:

        print("icType: ", icType)

        for res in reses:

            print("  Res: ", res)

            filenameOut = "particles_%s.nc" %res

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

                # create velocity and strain fields at material points
                # velocities
                uvVelMP = np.zeros([nParticles, 2])
                # strains
                strainAnalyticalMP = np.zeros([nParticles, 3])

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
                    u, v, strain11, strain22, strain12, divu, divv = velocities_strains_analytical(lat, lon, mu, lu, mv, lv)
                    uvVelMP[iParticle, 0] = u
                    uvVelMP[iParticle, 1] = v
                    strainAnalyticalMP[iParticle, 0] = strain11
                    strainAnalyticalMP[iParticle, 1] = strain22
                    strainAnalyticalMP[iParticle, 2] = strain12

                # append velocity and strain to particle file
                filenameOut = filenameIn
                fileOut = Dataset(filenameOut, "a", format="NETCDF3_CLASSIC")
                var = fileOut.createVariable("uvVelMP","d",dimensions=["nParticles", "TWO"])
                var[:] = uvVelMP[:]
                var = fileOut.createVariable("strainAnalyticalMP","d",dimensions=["nParticles", "THREE"])
                var[:] = strainAnalyticalMP[:]
                fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles()
