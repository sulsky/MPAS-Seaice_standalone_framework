import sys

sys.path.append("../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions

import os

#--------------------------------------------------------------------

def create_particles():

    reses = ["2562","10242","40962","163842"]

    icTypes = ["cosine_bell","slotted_cylinder"]

    for icType in icTypes:

        print("icType: ", icType)

        for res in reses:

            print("  Res: ", res)

            filenameOut = "particles_%s_%s.nc" %(icType,res)

            if (not os.path.isfile(filenameOut)):

                filenameMesh = "grid.%s.nc" %(res)
                filenameIceConc  = "ic_%s_%s.nc" %(icType,res)

                initial_particle_positions(filenameMesh,
                                           filenameIceConc,
                                           filenameOut,
                                           "number",
                                           9,
                                           "even",
                                           6371229.0)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles()
