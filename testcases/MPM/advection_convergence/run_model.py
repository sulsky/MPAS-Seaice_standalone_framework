import os
import time
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise
from datetime import datetime

#-------------------------------------------------------------------------------

def run_model(usePolympo):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    icTypes = ["cosine_bell","slotted_cylinder"]
    #icTypes = ["cosine_bell"]

    gridSizes = [2562, 10242, 40962, 163842]
    #gridSizes = [2562]

    for icType in icTypes:

         for gridSize in gridSizes:

             print()
             message = "IC type: %s, Gridsize: %i" %(icType, gridSize)
             print(message)
             print("-"*len(message))

             if (not os.path.isdir("output")):
                 os.mkdir("output")

             os.system("rm grid.nc ic.nc particles.nc")
             os.system("ln -s grid.%i.nc grid.nc" %(gridSize))
             os.system("ln -s ic_%s_%i.nc ic.nc" %(icType, gridSize))
             os.system("ln -s particles_%s_%i.nc particles.nc" %(icType, gridSize))

             os.system("rm -rf output_%s_%i" %(icType, gridSize))

             cmd = "%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE)
             start = time.perf_counter()
             os.system(cmd)
             end = time.perf_counter()
             print(f"Elapsed time: {end - start:.6f} seconds")

             os.system("mv output output_%s_%i" %(icType, gridSize))

             if (usePolympo):
                 usePolympoStr = "polympo"
             else:
                 usePolympoStr = "nonpolympo"
             os.system("mv log.seaice.0000.out log.seaice.0000.out_%s_%i_%s" %(icType, gridSize, usePolympoStr))

#-----------------------------------------------------------------------------

if __name__ == "__main__":

    run_model(False)
