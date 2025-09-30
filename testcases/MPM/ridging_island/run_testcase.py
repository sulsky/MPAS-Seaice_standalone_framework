import sys

sys.path.append("../../../utils/testcases/")
from create_square_quad_mesh import create_square_quad_mesh
from create_square_hex_mesh import create_square_hex_mesh
from plot_mesh import plot_mesh

sys.path.append("../../../utils/MPM/particle_initialization/")
from create_particles_from_cell_file import create_particles_from_cell_file

from create_ics import create_ics
from run_model import run_model
from plot_testcase import plot_testcase

import os
import numpy as np
import argparse
import subprocess
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def in_island(x,y):

    inIsland = False

    xy1 = 400000.0
    xy2 = 550000.0
    xy3 = 600000.0

    if ((x >= xy1 and x <= xy3 and y >= xy2 and y <= xy3) or
        (x >= xy2 and x <= xy3 and y >= xy1 and y <= xy3)):

        inIsland = True

    return inIsland

#-------------------------------------------------------------------------------

def cull_island(filenameIn,
                filenameOut):

    filein = Dataset(filenameIn,"a")
    nCells = len(filein.dimensions["nCells"])
    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]
    cullCell = filein.createVariable("cullCell","i",dimensions=["nCells"])

    for iCell in range(0,nCells):
        if (in_island(xCell[iCell],
                      yCell[iCell])):
            cullCell[iCell] = 1
        else:
            cullCell[iCell] = 0

    filein.close()

    MPAS_TOOLS_DIR = os.environ.get('MPAS_TOOLS_DIR')
    if (MPAS_TOOLS_DIR is None):
        raise Exception("MPAS_TOOLS_DIR environment variable must be set")

    MpasCellCuller = MPAS_TOOLS_DIR+"/mesh_tools/mesh_conversion_tools/MpasCellCuller.x"
    if (not os.path.isfile(MpasCellCuller)):
        raise Exception("MpasCellCuller executable must be built")

    subprocess.run([MpasCellCuller,filenameIn,filenameOut])

#-------------------------------------------------------------------------------

def run_testcase(runtype,
                 meshtype):

    print("Create grid...")
    lx = 1200000.0
    ly = 1200000.0

    x0 = 0.0
    y0 = 0.0

    if (meshtype == "quad"):

        nx = 120
        ny = 120

        gridFilenameNoIsland = create_square_quad_mesh(nx, ny,
                                                       lx, ly,
                                                       x0, y0)

    elif (meshtype == "hex"):

        dc = 10000.0

        gridFilenameNoIsland = create_square_hex_mesh(dc,
                                                      lx, ly,
                                                      x0, y0)

    else:
        raise Exception("Unknown mesh type: "+meshtype)

    # cull island
    gridFilename = os.path.splitext(os.path.basename(gridFilenameNoIsland))[0]+"_island.nc"
    cull_island(gridFilenameNoIsland,
                gridFilename)

    plot_mesh(gridFilename)

    print("Create ICs...")
    create_ics(gridFilename)

    print("Create particles...")
    particleInitType = "number"
    particleInitNumber = "9"
    particlePositionInitType = "even"
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

    parser = argparse.ArgumentParser()

    parser.add_argument('-r', dest="runtype",  choices=["orig","mpm"], default="orig")
    parser.add_argument('-m', dest="meshtype", choices=["quad","hex"], default="quad")

    args = parser.parse_args()

    run_testcase(args.runtype,
                 args.meshtype)
