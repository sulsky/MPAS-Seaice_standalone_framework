import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

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

            os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

            os.system("mv output output_%s_%i" %(test,gridSize))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
