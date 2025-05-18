import sys

sys.path.append("../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions

import os

#--------------------------------------------------------------------

def create_particles(res="2562"):

    #icTypes = ["cosine_bell","slotted_cylinder"]
    icTypes = ["slotted_cylinder"]

    for icType in icTypes:

        print("icType: ", icType)

        print("  Res: ", res)

        filenameOut = "particles_%s_%s.nc" %(icType,res)

        if (not os.path.isfile(filenameOut)):

            filenameMesh = "grid.%s.nc" %(res)

            initial_particle_positions(filenameMesh,
                                       filenameOut,
                                       "number",
                                       9,
                                       "onePerEdge",
                                       icType,
                                       6371229.0)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles(res)
