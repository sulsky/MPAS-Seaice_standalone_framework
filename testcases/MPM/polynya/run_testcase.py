import sys
sys.path.append("../../../utils/MPM/particle_initialization/")
from empty_particle_file import empty_particle_file

from create_ic import create_ic
from create_forcing import create_forcing
from plot_particles_scatter import plot_particles_scatter
from plot_cells_scatter import plot_cells_scatter
import os

#-------------------------------------------------------------------------------

def run_testcase():

    print("Get grid file")
    print("=============")
    MPAS_SEAICE_STANDALONE_DATA = os.environ.get('MPAS_SEAICE_STANDALONE_DATA')
    if (MPAS_SEAICE_STANDALONE_DATA is None):
        raise Exception("MPAS_SEAICE_STANDALONE_DATA must be set")
    filanameDst = "grid.nc"
    filenameSrc = "%s/testcases/polynya/grid_polynya.nc" %(MPAS_SEAICE_STANDALONE_DATA)
    if (not os.path.isfile(filanameDst)):
        os.symlink(filenameSrc, filanameDst)

    print("\nCreate empty particle file")
    print(  "==========================")
    empty_particle_file("particles.nc")

    print("\nCreate ICs")
    print(  "==========")
    create_ic()

    print("\nCreate forcing")
    print(  "==============")
    create_forcing()

    tests = ["cells",
             "particles"]

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    print("\nRun tests")
    print(  "=========")
    for test in tests:

        print("\nTest type: %s" %(test))
        print(  "-----------"+"-"*len(test))

        cmd = "rm namelist.seaice streams.seaice"
        print(cmd)
        os.system(cmd)

        cmd = "ln -s namelist.seaice.%s namelist.seaice" %(test)
        print(cmd)
        os.system(cmd)

        cmd = "ln -s streams.seaice.%s streams.seaice" %(test)
        print(cmd)
        os.system(cmd)

        cmd = "rm log.seaice.0000.out"
        print(cmd)
        os.system(cmd)

        cmd = "rm -rf output_%s" %(test)
        print(cmd)
        os.system(cmd)

        cmd = "%s" %(MPAS_SEAICE_EXECUTABLE)
        print(cmd)
        os.system(cmd)

        cmd = "cp log.seaice.0000.out log.seaice.0000.out_%s" %(test)
        print(cmd)
        os.system(cmd)

        if (test == "cells"):
            plot_cells_scatter()
        elif (test == "particles"):
            plot_particles_scatter()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
