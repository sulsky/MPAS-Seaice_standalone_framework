import os
import argparse

#-------------------------------------------------------------------------------

def run_model(nProcs):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    if (not os.path.isdir("output")):
        os.mkdir("output")

    os.system("rm -rf output_%i" %(nProcs))

    os.system("mpirun -np %i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE))

    os.system("mv output output_%i" %(nProcs))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', dest="nProcs", type=int)
    args = parser.parse_args()

    run_model(args.nProcs)
