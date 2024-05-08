import sys

sys.path.append("../../advection")

from get_testcase_data import get_testcase_data
from create_ics import create_ics
from create_particles import create_particles
from run_model import run_model
from plot_testcase import plot_testcase
from advection_map import advection_map
from advection_equatorial import advection_equatorial
from advection_error_convergence import advection_error_convergence

get_testcase_data()

create_ics()

create_particles()

#run_model()

plot_testcase()

advection_map()

advection_equatorial()

advection_error_convergence()
