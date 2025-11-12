from run_model import run_model

import sys, os

sys.path.append("../../../../utils/MPM/particle_initialization")
from empty_particle_file import empty_particle_file

sys.path.append("../../../../utils/testcases")
from create_square_hex_mesh import create_square_hex_mesh
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
    lx = 1.0
    ly = 1.0
    x0 = 0.0
    y0 = 0.0
    dc = 0.05
    gridFilename = create_square_hex_mesh(dc,lx,ly,x0,y0)
    plot_mesh(gridFilename)
    print("Sym link grid...")
    if (os.path.isfile("grid.nc")):
        os.system("rm grid.nc")
    os.symlink(gridFilename,"grid.nc")

    print("Create particle file")
    print("=================")
    empty_particle_file('particles.nc')

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
