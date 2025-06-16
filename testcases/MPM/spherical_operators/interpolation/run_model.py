import sys

sys.path.append("../../../utils/testcases")
from log_messages import log_message
from execute_model import execute_model

import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model(logFile):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    gridSizes = [2562, 10242, 40962, 163842]

    tests = ["1-x", "y-z", "latlon", "nonlin"]
    #tests = ["1-x"]

    for test in tests:

        print(" test: ", test)

        for gridSize in gridSizes:

            print("  Gridsize: ", gridSize)

            os.system("rm grid.nc ic.nc particles.nc")
            os.system("ln -s grid.%i.nc grid.nc" %(gridSize))
            os.system("ln -s ic_%s_%i.nc ic.nc" %(test,gridSize))
            os.system("ln -s particles_%s_%i.nc particles.nc" %(test,gridSize))

            os.system("rm -rf output_%s_%i" %(test,gridSize))

            execute_model(MPAS_SEAICE_EXECUTABLE,
                          1,
                          logFile)

            os.system("mv output output_%s_%i" %(test,gridSize))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
