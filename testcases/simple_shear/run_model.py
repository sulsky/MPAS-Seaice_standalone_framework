import os

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    if (not os.path.isdir("output")):
        os.mkdir("output")

    os.system("ln -s namelist.seaice.simple_shear namelist.seaice")
    os.system("ln -s streams.seaice.simple_shear streams.seaice")

    os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
