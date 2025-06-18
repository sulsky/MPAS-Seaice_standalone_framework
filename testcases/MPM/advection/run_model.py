import sys

sys.path.append("../../../utils/testcases")
from execute_model import execute_model

import os
import argparse

#-------------------------------------------------------------------------------

def run_model(nCells,
              nProcs,
              runtype,
              logFile):

    print()
    message = "Run model for nCells: %i and nProcs: %i" %(nCells, nProcs)
    print(message)
    print("-"*len(message))

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    cmd = "rm grid.nc ic.nc particles.nc log.seaice.*"
    print(cmd)
    os.system(cmd)

    cmd = "rm -rf output output_%s_%i" %(runtype,nProcs)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s grid.%i.nc grid.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s ic_slotted_cylinder_%i.nc ic.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s particles_slotted_cylinder_%i.nc particles.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    if (nProcs > 1):
        os.chdir("./graphs")

        cmd = "ln -s graph.%i.info.part.%i graph.info.part.%i" %(nCells, nProcs, nProcs)
        print(cmd)
        os.system(cmd)

        os.chdir("..")

    execute_model(MPAS_SEAICE_EXECUTABLE,
                  nProcs,
                  logFile)

    cmd = "mv output output_%s_%i" %(runtype, nProcs)
    print(cmd)
    os.system(cmd)

    print()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest="nCells", type=int, required=True)
    parser.add_argument('-p', dest="nProcs", type=int, required=True)
    parser.add_argument('-r', dest="runtype", type=char, required=False)
    args = parser.parse_args()

    run_model(args.nCells, args.nProcs, args.runtype)
