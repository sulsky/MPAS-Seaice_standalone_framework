import os, sys

sys.path.append("../../../../utils/testcases")
from execute_model import execute_model

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    if (os.path.isdir("output")):
        cmd = "rm -rf output"
        os.system(cmd)
    os.mkdir("output")

    nProcs = 1
    logFile = None
    execute_model(MPAS_SEAICE_EXECUTABLE,
                  nProcs,
                  logFile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
