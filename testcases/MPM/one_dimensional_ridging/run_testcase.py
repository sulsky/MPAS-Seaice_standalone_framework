import sys

sys.path.append("../../../utils/testcases/")
from create_square_quad_mesh import create_square_quad_mesh

sys.path.append("../../../utils/MPM/particle_initialization/")
from create_particles_from_cell_file import create_particles_from_cell_file

from create_ics import create_ics
from run_model import run_model
from plot_heatmap import plot_heatmap

import os

#-------------------------------------------------------------------------------

def run_testcase():

    runtype = "orig"

    nx = 100
    ny = 100

    lx = 1000000.0
    ly = 1000000.0

    x0 = 0.0
    y0 = 0.0

    particleInitType = "number"
    particleInitNumber = "9"
    particlePositionInitType = "even"

    print("Create grid...")
    gridFilename = create_square_quad_mesh(nx, ny,
                                           lx, ly,
                                           x0, y0)

    print("Create ICs...")
    create_ics(gridFilename)

    print("Create particles...")
    create_particles_from_cell_file(gridFilename,
                                    "ic.nc",
                                    particleInitType,
                                    particleInitNumber,
                                    particlePositionInitType,
                                    "particles.nc")

    print("Sym link grid...")
    if (os.path.isfile("grid.nc")):
        os.system("rm grid.nc")
    os.symlink(gridFilename,"grid.nc")

    print("Run model...")
    run_model(runtype)

    print("Plot heatmap...")
    plot_heatmap("grid.nc",
                 "output_%s/output.2000.nc" %(runtype))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
