import sys

sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

sys.path.append("../../advection")
from create_ics import create_ics

sys.path.append("../../../utils/MPM/particle_initialization/")
from empty_particle_file import empty_particle_file
from run_model import run_model
from check_results import check_results

nCells = 10242
n1 = 16
n2 = 32

get_testcase_data_spherical()

create_ics(nCells)

empty_particle_file("particles.nc")

run_model(nCells, n1)

run_model(nCells, n2)

check_results(n1, n2)
