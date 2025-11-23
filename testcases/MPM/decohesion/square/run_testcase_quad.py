from run_model import run_model
from create_ics import create_ics

import sys, os

sys.path.append("../../../../utils/MPM/particle_initialization/")
from create_particles_from_cell_file import create_particles_from_cell_file
from initial_particle_positions import initial_particle_positions

sys.path.append("../../../../utils/testcases")
from create_square_quad_mesh import create_square_quad_mesh
from plot_mesh import plot_mesh

import argparse

#-------------------------------------------------------------------------------

def run_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nDecohesion test case\n")
        logFile.write(  "===================\n")
        logFile.flush()
    else:
        logFile = None

    print("Create grid")
    print("=================")
    lx = 240000.0
    ly = 240000.0
    x0 = 0.0
    y0 = 0.0
    nx = 24
    ny = 24

    gridFilename = create_square_quad_mesh(nx, ny,
                                           lx, ly,
                                           x0, y0)

    plot_mesh(gridFilename)

    print("Sym link grid...")
    if (os.path.isfile("grid.nc")):
        os.system("rm grid.nc")
    os.symlink(gridFilename,"grid.nc")

    print("Create ICs...")
    create_ics(gridFilename)

    print("Create particles...")
    particleInitType = "number"
    particleInitNumber = "4"
    particlePositionInitType = "even"
    particleGeometry = "square"
    sphereRadius = 1.0
    initial_particle_positions(gridFilename,
                                "particles.nc",
                                particleInitType,
                                particleInitNumber,
                                particlePositionInitType,
                                particleGeometry,
                                sphereRadius)

    print("Run model")
    print("=================")
    run_model()

    if (logFile is not None):
        logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.logFilename)
