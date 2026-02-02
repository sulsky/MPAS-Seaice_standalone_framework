from run_model import run_model
from create_ics import create_ics
from plot_decohesion import plot_decohesion

import sys, os

sys.path.append("../../../../utils/MPM/particle_initialization/")
from initial_particle_positions import initial_particle_positions

sys.path.append("../../../../utils/testcases")
from create_square_hex_mesh import create_square_hex_mesh
from create_square_quad_mesh import create_square_quad_mesh
from plot_mesh import plot_mesh

import argparse

#-------------------------------------------------------------------------------

def run_testcase(meshtype, testType, logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nDecohesion test case\n")
        logFile.write(  "===================\n")
        logFile.flush()
    else:
        logFile = None

    print("Create grid")
    print("=================")
    if (meshtype == 'hex'):
        lx = 240000.0
        ly = 240000.0
        x0 = 0.0
        y0 = 0.0
        dc = 10000.0
        gridFilename = create_square_hex_mesh(dc,
                                          lx, ly,
                                          x0, y0)
    elif (meshtype == 'quad'):
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
    print("=================")
    create_ics(gridFilename, testType)

    print("Create particles...")
    print("=================")
    if (meshtype == 'hex'):
       particleInitType = "number"
       particleInitNumber = "9"
       particlePositionInitType = "onePerEdge"
       particleGeometry = "square"
       sphereRadius = 1.0
       initial_particle_positions(gridFilename,
                                "particles.nc",
                                particleInitType,
                                particleInitNumber,
                                particlePositionInitType,
                                particleGeometry,
                                sphereRadius)
    elif (meshtype == 'quad'):
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

    print("Plot decohesion")
    print("=================")
    plot_decohesion(gridFilename)

    if (logFile is not None):
        logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest="meshtype", choices=["quad","hex"], default="quad")
    parser.add_argument('-t', dest="testType", choices=["tension","shear"], default="tension")
    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.meshtype, args.testType, args.logFilename)
