import sys

sys.path.append("../../../utils/testcases")
from execute_model import execute_model

import os
import argparse

#-------------------------------------------------------------------------------

def run_model(runtype,
              nProcs=1):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    cmd = "rm namelist.seaice streams.seaice"
    print(cmd)
    os.system(cmd)

    cmd = "ln -s namelist.seaice.%s namelist.seaice" %(runtype)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s streams.seaice.%s streams.seaice" %(runtype)
    print(cmd)
    os.system(cmd)

    cmd = "rm -rf output output_%s restarts" %(runtype)
    print(cmd)
    os.system(cmd)

    execute_model(MPAS_SEAICE_EXECUTABLE,
                  nProcs)

    cmd = "mv output output_%s" %(runtype)
    print(cmd)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-r', dest="runtype", choices=["orig","mpm"], required=True)
    parser.add_argument('-n', dest="nprocs", type=int, default=1)

    args = parser.parse_args()

    run_model(args.runtype,
              args.nprocs)
