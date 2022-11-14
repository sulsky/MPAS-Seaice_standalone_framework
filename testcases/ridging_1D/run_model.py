import os

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    if (not os.path.isdir("output")):
        os.mkdir("output")

    os.system("rm -rf namelist.seaice streams.seaice")
    os.system("ln -s namelist.seaice.ridging_1D namelist.seaice")
    os.system("ln -s streams.seaice.ridging_1D streams.seaice")

    os.system("%s" %(MPAS_SEAICE_EXECUTABLE))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
