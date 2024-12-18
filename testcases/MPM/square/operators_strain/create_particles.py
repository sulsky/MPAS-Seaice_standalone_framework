from netCDF4 import Dataset
import numpy as np
from math import sin, cos, pi
import sys
sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions
sys.path.append("../../../square/operators_strain")
from create_ics import velocities_strains_stress_divergences

#-------------------------------------------------------------

def particles_ic(gridfile, icfile, velType):

    # load grid file
    grid = Dataset(gridfile, "r")

    nVertices = len(grid.dimensions["nVertices"])
    xVertex = grid.variables["xVertex"][:]
    yVertex = grid.variables["yVertex"][:]

    grid.close()

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    # create material points
    initial_particle_positions(gridfile,
                               icfile,
                               "number",
                               9,
                               "even",
                               "uniform",
                               1.0)

    # read back in particle positions
    filenameIn = icfile
    fileIn = Dataset(filenameIn, "r")
    nParticles = len(fileIn.dimensions["nParticles"])
    posnMP = fileIn.variables["posnMP"][:]
    fileIn.close()

    # calculate output variables
    uvVelMP = np.zeros([nParticles, 2])
    # strains
    strainAnalyticalMP = np.zeros([nParticles, 3])
    # stress
    stressAnalyticalMP = np.zeros([nParticles, 3])

    for iParticle in range(0, nParticles):
        xp = posnMP[iParticle, 0]
        yp = posnMP[iParticle, 1]
        zp = posnMP[iParticle, 2]

        x = xp - xMin
        y = yp - yMin

        u, v, e11, e22, e12, divu, divv = velocities_strains_stress_divergences(x, y, velType)

        uvVelMP[iParticle, 0] = u
        uvVelMP[iParticle, 1] = v
        strainAnalyticalMP[iParticle, 0] = e11
        strainAnalyticalMP[iParticle, 1] = e22
        strainAnalyticalMP[iParticle, 2] = e12

        # Linear constitutive model
        lam = 1.0
        stressAnalyticalMP[iParticle, 0] = lam*e11
        stressAnalyticalMP[iParticle, 1] = lam*e22
        stressAnalyticalMP[iParticle, 2] = lam*e12

    # append velocity, stress, strain to particle file
    fileOut = Dataset(icfile, "a", format="NETCDF3_CLASSIC")
    var = fileOut.createVariable("uvVelMP","d",dimensions=["nParticles", "TWO"])
    var[:] = uvVelMP[:]
    var = fileOut.createVariable("strainAnalyticalMP","d",dimensions=["nParticles", "THREE"])
    var[:] = strainAnalyticalMP[:]
    var = fileOut.createVariable("stressAnalyticalMP","d",dimensions=["nParticles", "THREE"])
    var[:] = stressAnalyticalMP[:]
    fileOut.close()

#-------------------------------------------------------------

def create_particles():

    gridTypes = ["hex","quad"]

    grids = {"hex": ["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}

    velType = "sinusoid"

    for gridType in gridTypes:
        for grid in grids[gridType]:

            gridfile = "grid_%s_%s.nc" %(gridType,grid)
            icfile = "particles_%s_%s.nc" %(gridType,grid)

            particles_ic(gridfile, icfile, velType)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles()
