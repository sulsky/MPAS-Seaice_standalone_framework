from run_model import run_model
from create_ics import create_ics
from add_initial_area_volume_categories import add_initial_area_volume_categories
from plot_results import plot_results

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
        logFile.write("\nDecohesion test case\n")
        logFile.write(  "===================\n")
        logFile.flush()
    else:
        logFile = None

    print("Create grid")
    print("=================")
    if (meshtype == 'hex'):
        lx = 1040000.0
        ly = 1040000.0
        x0 = 0.0
        y0 = 0.0
        dc = 10400.0
        gridFilename = create_square_hex_mesh(dc,
                                          lx, ly,
                                          x0, y0)
    elif (meshtype == 'quad'):
        lx = 1040000.0
        ly = 1040000.0
        x0 = 0.0
        y0 = 0.0
        nx = 104
        ny = 104
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
    create_ics(gridFilename)

    print("Create particles...")
    print("=================")
    if (meshtype == 'hex'):
       particleInitType = "number"
       particleInitNumber = "9"
       particlePositionInitType = "onePerEdge"
       particleGeometry = "BBMsquare"
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
       particleGeometry = "BBMsquare"
       sphereRadius = 1.0
       initial_particle_positions(gridFilename,
                                "particles.nc",
                                particleInitType,
                                particleInitNumber,
                                particlePositionInitType,
                                particleGeometry,
                                sphereRadius)

    print("Adding ice categories...")
    print("=================")
    add_initial_area_volume_categories()

    print("Run model...")
    print("=================")
    run_model()

    print("Plot results...")
    print("=================")
    plot_results(gridFilename)

    if (logFile is not None):
        logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest="meshtype", choices=["quad","hex"], default="hex")
    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.meshtype, args.logFilename)
