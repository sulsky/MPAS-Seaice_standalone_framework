from run_model import run_model
from create_particles import create_particles

import sys, os

sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions

sys.path.append("../../../../utils/testcases")
from create_square_hex_mesh import create_square_hex_mesh
from create_square_quad_mesh import create_square_quad_mesh
from plot_mesh import plot_mesh

import argparse

#-------------------------------------------------------------------------------

def run_testcase(meshtype, logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nElastic Disks test case\n")
        logFile.write(  "===================\n")
        logFile.flush()
    else:
        logFile = None

    print("Create grid")
    print("=================")
    if (meshtype == 'hex'):
        lx = 1.0
        ly = 1.0
        x0 = 0.0
        y0 = 0.0
        dc = 0.05
        gridFilename = create_square_hex_mesh(dc,
                                          lx, ly,
                                          x0, y0)
    elif (meshtype == 'quad'):
        lx = 1.0
        ly = 1.0
        x0 = 0.0
        y0 = 0.0
        nx = 20
        ny = 20
        gridFilename = create_square_quad_mesh(nx, ny,
                                           lx, ly,
                                           x0, y0)

    plot_mesh(gridFilename)

    print("Sym link grid...")
    if (os.path.isfile("grid.nc")):
        os.system("rm grid.nc")
    os.symlink(gridFilename,"grid.nc")

    print("Create particles...")
    print("=================")
    create_particles(gridFilename, meshtype)

    print("Run model")
    print("=================")
    run_model()

    if (logFile is not None):
        logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest="meshtype", choices=["quad","hex"], default="quad")
    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.meshtype, args.logFilename)
