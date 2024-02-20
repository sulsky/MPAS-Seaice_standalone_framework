import sys

sys.path.append("../../advection")
from get_testcase_data import get_testcase_data
from create_ics import create_ics
from create_particles import create_particles

sys.path.append("./graphs")
from create_graph_file_basic import create_graph_file_basic
from create_graph_file_metis import create_graph_file_metis

from check_particle_positions_start_end import check_particle_positions_start_end
from check_particle_positions_nprocs import check_particle_positions_nprocs
from run_model import run_model
import os

nCells = 2562

get_testcase_data()

create_ics()

create_particles()

create_graph_file_basic(nCells, 2)
create_graph_file_basic(nCells, 4)

create_graph_file_metis(nCells, 16)
create_graph_file_metis(nCells, 32)

run_model(nCells, 1)

check_particle_positions_start_end("./output_1/particles*")

run_model(nCells, 2)
check_particle_positions_nprocs(2,1)

run_model(nCells, 4)
check_particle_positions_nprocs(4,1)

run_model(nCells, 16)
check_particle_positions_nprocs(16,1)

run_model(nCells, 32)
check_particle_positions_nprocs(32,1)
check_particle_positions_nprocs(32,16)
