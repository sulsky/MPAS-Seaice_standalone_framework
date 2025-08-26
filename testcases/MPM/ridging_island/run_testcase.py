import sys

sys.path.append("../../../utils/testcases/")
from create_square_quad_mesh import create_square_quad_mesh

sys.path.append("../../../utils/MPM/particle_initialization/")
from create_particles_from_cell_file import create_particles_from_cell_file

from create_ics import create_ics
from run_model import run_model
from plot_testcase import plot_testcase

import os
import numpy as np

#-------------------------------------------------------------------------------

def run_testcase():

    runtype = "orig"

    nx = 120
    ny = 120

    lx = 1200000.0
    ly = 1200000.0

    x0 = 0.0
    y0 = 0.0

    particleInitType = "number"
    particleInitNumber = "9"
    particlePositionInitType = "even"

    # island
    cull = np.zeros((nx,ny),dtype="i")
    for ix in range(0,nx):
        for iy in range(0,ny):
            if ((ix >= 40 and ix < 60 and iy >= 55 and iy < 60) or
                (ix >= 55 and ix < 60 and iy >= 40 and iy < 60)):
                cull[ix,iy] = 1

    print("Create grid...")
    gridFilename = create_square_quad_mesh(nx, ny,
                                           lx, ly,
                                           x0, y0,
                                           cull)

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

    # plot output
    plot_testcase("./output_%s/output.2000.nc" %(runtype))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
