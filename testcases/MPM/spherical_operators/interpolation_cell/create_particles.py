from netCDF4 import Dataset
import numpy as np
import math
import sys
from create_ic import latlon_from_xyz, grid_rotation_forward, icearea_analytical
sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions
import os

#--------------------------------------------------------------------

def create_particles():

    rotateCartesianGrid = True
    r = 1.0

    reses = ["2562","10242","40962","163842"]

    icTypes = ["uniform"]

    tests = ["1", "1+x", "1+y", "1+z", "lat", "lon", "nonlin"]

    rotateCartesianGrid = True
    r = 1.0

    for test in tests:
        if (test == "1"):
           c = [1, 0, 0, 0, 0, 0, 0]
        elif (test == "1+x"):
           c = [1, 1, 0, 0, 0, 0, 0]
        elif (test == "1+y"):
           c = [1, 0, 1, 0, 0, 0, 0]
        elif (test == "1+z"):
           c = [1, 0, 0, 1, 0, 0, 0]
        elif (test == "lat"):
           c = [0, 0, 0, 0, 1, 0, 0]
        elif (test == "lon"):
           c = [0, 0, 0, 0, 0, 1, 0]
        elif (test == "nonlin"):
           c = [0, 0, 0, 0, 0, 0, 1]

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

                    # read position data and iceAreaCell from particle file
                    fileOut = Dataset(filenameOut, "r+", format="NETCDF3_CLASSIC")
                    nParticles = len(fileOut.dimensions["nParticles"])
                    posnMP = fileOut.variables["posnMP"][:]
                    iceAreaCellMP = fileOut.variables["iceAreaCellMP"][:]
                    iceAreaCategoryMP = fileOut.variables["iceAreaCategoryMP"][:]

                    # modify iceAreaCell at material points
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
                        icearea = icearea_analytical(xp, yp, zp, lat, lon, c)
                        iceAreaCellMP[iParticle] = icearea
                        iceAreaCategoryMP[iParticle, 0] = icearea

                    fileOut.variables["iceAreaCellMP"][:] = iceAreaCellMP
                    fileOut.variables["iceAreaCategoryMP"][:] = iceAreaCategoryMP

                    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles()
