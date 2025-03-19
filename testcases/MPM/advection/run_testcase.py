import sys

sys.path.append("../../advection")
from get_testcase_data import get_testcase_data
from create_ics import create_ics
from add_deldyn_to_ics import add_deldyn_to_ics
from create_particles import create_particles

sys.path.append("./graphs")
from create_graph_file_basic import create_graph_file_basic
from create_graph_file_metis import create_graph_file_metis

from check_particle_positions_start_end import check_particle_positions_start_end
from check_particle_positions_nprocs import check_particle_positions_nprocs
from check_particles_moved import check_particles_moved
from run_model import run_model
import os

nCells = 2562

particleFilename1 = "particles_output.2000-01-01_00.00.00.nc"
particleFilename2 = "particles_output.2000-03-01_00.00.00.nc"

get_testcase_data()

create_ics()

add_deldyn_to_ics(3600.0)

create_particles()

create_graph_file_basic(nCells, 2)
create_graph_file_basic(nCells, 4)

create_graph_file_metis(nCells, 16)
create_graph_file_metis(nCells, 32)

run_model(nCells, 1)
check_particles_moved("./output_1/"+particleFilename1,
                      "./output_1/"+particleFilename2)

check_particle_positions_start_end("./output_1/particles*")

run_model(nCells, 2)
check_particles_moved("./output_2/"+particleFilename1,
                      "./output_2/"+particleFilename2)
check_particle_positions_nprocs(2,1)

run_model(nCells, 4)
check_particles_moved("./output_4/"+particleFilename1,
                      "./output_4/"+particleFilename2)
check_particle_positions_nprocs(4,1)

run_model(nCells, 16)
check_particles_moved("./output_16/"+particleFilename1,
                      "./output_16/"+particleFilename2)
check_particle_positions_nprocs(16,1)

run_model(nCells, 32)
check_particles_moved("./output_32/"+particleFilename1,
                      "./output_32/"+particleFilename2)
check_particle_positions_nprocs(32,1)
check_particle_positions_nprocs(32,16)
