import sys

sys.path.append("../../../utils/testcases/")
from create_square_quad_mesh import create_square_quad_mesh

sys.path.append("../../../utils/MPM/particle_initialization/")
from create_particles_from_cell_file import create_particles_from_cell_file

from create_forcing import create_forcing
from create_ics import create_ics
from plot_heatmap import plot_heatmap

import os

#-------------------------------------------------------------------------------

def run_testcase():

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

    print("Create forcing...")
    create_forcing(gridFilename)

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
    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))
    os.system(MPAS_SEAICE_EXECUTABLE)

    print("Plot heatmap...")
    plot_heatmap("grid.nc",
                 "output/output.2000.nc")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
