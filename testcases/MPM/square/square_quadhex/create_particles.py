import sys
import os

sys.path.append("../../../../utils/MPM/particle_initialization")
from create_particles_from_cell_file import create_particles_from_cell_file

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

    for gridType in gridTypes:
        for grid in grids[gridType]:

            gridfile = "grid_%s_%s.nc" %(gridType,grid)
            icfile = "ic_%s_%s.nc" %(gridType,grid)
            partfile = "particles_%s_%s.nc" %(gridType,grid)

            if not os.path.isfile(partfile):
               create_particles_from_cell_file(gridfile,
                                               icfile,
                                               'number',
                                               9,
                                               'even',
                                               partfile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_particles()
