import sys

sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

sys.path.append("../../advection")
from create_ics import create_ic_file
from add_deldyn_to_ics import add_deldyn_to_ics
from create_particles import create_particles

sys.path.append("./graphs")
from create_graph_file_basic import create_graph_file_basic
from create_graph_file_metis import create_graph_file_metis

sys.path.append("../../../testing")
from testing_utils import create_new_namelist

from check_particle_positions_start_end import check_particle_positions_start_end
from check_particle_positions_nprocs import check_particle_positions_nprocs
from check_particles_moved import check_particles_moved
from run_model import run_model

import math
import argparse
import os

#-------------------------------------------------------------------------------

def run_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nAdvection test case\n")
        logFile.write(  "===================\n")
        logFile.flush()

    nCells = 2562

    particleFilename1 = "particles_output.2000-01-01_00.00.00.nc"
    particleFilename2 = "particles_output.2000-01-11_00.00.00.nc"

    runtypes = [{"name":"original",
                 "polympo":False,
                 "nProcs":[{"nProcsRun":2,  "nProcsComp":[1]},
                           {"nProcsRun":4,  "nProcsComp":[1]},
                           {"nProcsRun":16, "nProcsComp":[1]},
                           {"nProcsRun":32, "nProcsComp":[1,16]}]},
                {"name":"polympo",
                 "polympo":True,
                 "nProcs":[]}]

    print("Get testcase data")
    print("=================")
    get_testcase_data_spherical()

    print("\nCreate ICs")
    print("==========")
    create_ic_file("2562", "slotted_cylinder" , math.pi / 6.0)
    create_ic_file("2562", "cosine_bell"      , math.pi / 6.0)

    print("\nAdd deldyn to ICs")
    print("=================")
    add_deldyn_to_ics(3600.0)

    print("\nCreate particles")
    print("================")
    create_particles()

    print("\nCreate graph files")
    print("==================")
    create_graph_file_basic(nCells, 2)
    create_graph_file_basic(nCells, 4)

    create_graph_file_metis(nCells, 16)
    create_graph_file_metis(nCells, 32)

    print("\nRun models")
    print("==========")
    for runtype in runtypes:
        print("Runtype name: ",runtype["name"])

        nmlChanges = {"mpm":{"config_use_mpm_polympo":runtype["polympo"]}}
        create_new_namelist("namelist.seaice.advection", "namelist.seaice", nmlChanges)

        run_model(nCells, 1, runtype["name"])
        check_particles_moved(("./output_%s_1/"%(runtype["name"]))+particleFilename1,
                              ("./output_%s_1/"%(runtype["name"]))+particleFilename2,
                              logFile)

        check_particle_positions_start_end("./output_%s_1/particles*"%(runtype["name"]),
                                           logFile)

        for nProcs in runtype["nProcs"]:

            print("nProcsRun:",nProcs["nProcsRun"])

            run_model(nCells, nProcs["nProcsRun"], runtype["name"])
            check_particles_moved(("./output_%s_%i/" %(runtype["name"],nProcs["nProcsRun"]))+particleFilename1,
                                  ("./output_%s_%i/" %(runtype["name"],nProcs["nProcsRun"]))+particleFilename2,
                                  logFile)

            for nProcsComp in nProcs["nProcsComp"]:
                check_particle_positions_nprocs(nProcs["nProcsRun"],
                                                nProcsComp,
                                                runtype["name"],
                                                logFile)

    if (logFilename is not None):
        logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.logFilename)
