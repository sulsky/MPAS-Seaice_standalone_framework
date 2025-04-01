import sys

sys.path.append("../../advection")
from get_testcase_data import get_testcase_data
from create_ics import create_ics
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
import os

nCells = 2562

particleFilename1 = "particles_output.2000-01-01_00.00.00.nc"
particleFilename2 = "particles_output.2000-03-01_00.00.00.nc"

runtypes = [{"name":"original",
             "polympo":False,
             "nProcs":[{"nProcsRun":2,  "nProcsComp":[1]},
                       {"nProcsRun":4,  "nProcsComp":[1]},
                       {"nProcsRun":16, "nProcsComp":[1]},
                       {"nProcsRun":32, "nProcsComp":[1,16]}]},
            {"name":"polympo",
             "polympo":True,
             "nProcs":[]}]

get_testcase_data()

create_ics()

add_deldyn_to_ics(3600.0)

create_particles()

create_graph_file_basic(nCells, 2)
create_graph_file_basic(nCells, 4)

create_graph_file_metis(nCells, 16)
create_graph_file_metis(nCells, 32)

for runtype in runtypes:
    print("Runtype name: ",runtype["name"])

    nmlChanges = {"mpm":{"config_use_mpm_polympo":runtype["polympo"]}}
    create_new_namelist("namelist.seaice.advection", "namelist.seaice", nmlChanges)

    run_model(nCells, 1)
    check_particles_moved("./output_1/"+particleFilename1,
                          "./output_1/"+particleFilename2)

    check_particle_positions_start_end("./output_1/particles*")

    for nProcs in runtype["nProcs"]:

        print("nProcsRun:",nProcs["nProcsRun"])

        run_model(nCells, nProcs["nProcsRun"])
        check_particles_moved(("./output_%i/" %(nProcs["nProcsRun"]))+particleFilename1,
                              ("./output_%i/" %(nProcs["nProcsRun"]))+particleFilename2)

        for nProcsComp in nProcs["nProcsComp"]:
            check_particle_positions_nprocs(nProcs["nProcsRun"],nProcsComp)
