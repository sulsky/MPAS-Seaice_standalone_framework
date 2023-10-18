import sys
sys.path.append("../../advection")
from get_testcase_data import get_testcase_data
from create_ics import create_ics
from check_particle_positions import check_particle_positions
from run_model import run_model
import os

nCells = 10242

get_testcase_data()

create_ics()

run_model(nCells)

check_particle_positions()
