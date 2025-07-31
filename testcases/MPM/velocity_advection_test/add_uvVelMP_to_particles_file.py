import argparse
from netCDF4 import Dataset
import math
import os

from create_ics import latlon_vector_rotation_forward, latlon_from_xyz

#--------------------------------------------------------------------

def add_uvVelMP_to_particles_file(earthRadius, rotateCartesianGrid):

    reses = ["2562","10242","40962","163842"]

    for res in reses:

       print("  Res: ", res)

       icFilename = "particles_%s.nc" %(res)

       if (not os.path.isfile(icFilename)):
           raise Exception("IC file missing: "+icFilename)

       filein = Dataset(icFilename,"a")

       try:
           filein.createDimension("TWO",2)
       except:
           pass

       try:
          uvVelMP = filein.createVariable("uvVelMP","d",dimensions=["nParticles","TWO"])
       except:
          uvVelMP = filein.variables["uvVelMP"][:]

       try:
          latCellMP = filein.createVariable("latCellMP","d",dimensions=["nParticles"])
       except:
          latCellMP = filein.variables["latCellMP"][:]

       try:
          lonCellMP = filein.createVariable("lonCellMP","d",dimensions=["nParticles"])
       except:
          lonCellMP = filein.variables["lonCellMP"][:]

       nParticles = len(filein.dimensions["nParticles"])
       posnMP = filein.variables["posnMP"][:,:]

       for iParticle in range(0, nParticles):

            latCellMP[iParticle], lonCellMP[iParticle] = latlon_from_xyz(
                                                           posnMP[iParticle, 0],
                                                           posnMP[iParticle, 1],
                                                           posnMP[iParticle, 2],
                                                           earthRadius)

            uvVelMP[iParticle, 0] = 40.0 * math.cos(latCellMP[iParticle])
            uvVelMP[iParticle, 1] = 0.0

            uvVelMP[iParticle, 0], uvVelMP[iParticle, 1] = latlon_vector_rotation_forward(
                                                           uvVelMP[iParticle, 0],
                                                           uvVelMP[iParticle, 1],
                                                           latCellMP[iParticle],
                                                           lonCellMP[iParticle],
                                                           posnMP[iParticle, 0],
                                                           posnMP[iParticle, 1],
                                                           posnMP[iParticle, 2],
                                                           earthRadius,
                                                           rotateCartesianGrid)

       filein.close()

#--------------------------------------------------------------------

if __name__ == "__main__":

    add_uvVelMP_to_particles_file(earthRadius, rotateCartesianGrid)
