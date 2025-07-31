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
    message = "Run model for nCells: %s and nProcs: %i" %(nCells, nProcs)
    print(message)
    print("-"*len(message))

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    icTypes = ["cosine_bell","slotted_cylinder"]
    #icTypes = ["cosine_bell"]

    tsteps = [3600, 1800, 900, 450]

    for icType in icTypes:

        for tstep in tsteps:

            print()
            message = "IC type: %s, gridSize: %s, Time Step: %i" %(icType, nCells, tstep)
            print(message)
            print("-"*len(message))

            if (not os.path.isdir("output")):
                os.mkdir("output")

            os.system("rm grid.nc ic.nc particles.nc namelist.seaice")
            os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(runtype, tstep))
            os.system("ln -s grid.%s.nc grid.nc" %(nCells))
            os.system("ln -s ic_%s_%i.nc ic.nc" %(icType, tstep))
            os.system("ln -s particles_%s_%s.nc particles.nc" %(icType, nCells))

            os.system("rm -rf output_%s_%i_%s" %(icType, tstep, runtype))


            if (nProcs > 1):
                os.chdir("./graphs")

                cmd = "ln -s graph.%s.info.part.%i graph.info.part.%i" %(nCells, nProcs, nProcs)
                print(cmd)
                os.system(cmd)

                os.chdir("..")

            execute_model(MPAS_SEAICE_EXECUTABLE,
                          nProcs,
                          logFile)

            os.system("mv output output_%s_%i_%s" %(icType, tstep, runtype))
            os.system("mv log.seaice.0000.out log.seaice.0000.out_%s_%i_%s" %(icType, tstep, runtype))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest="nCells", type=char, required=True)
    parser.add_argument('-p', dest="nProcs", type=int, required=True)
    parser.add_argument('-r', dest="runtype", type=char, required=False)
    args = parser.parse_args()

    run_model(args.nCells, args.nProcs, args.runtype)
