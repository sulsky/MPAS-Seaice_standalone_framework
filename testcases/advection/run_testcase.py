import sys

sys.path.append("../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from create_ics import create_ics
from run_model import run_model
from aggregate_tracers import aggregate_tracers
from plot_testcase import plot_testcase
from advection_map import advection_map
from advection_equatorial import advection_equatorial
from advection_error_convergence import advection_error_convergence

get_testcase_data_spherical()

create_ics()

run_model()

aggregate_tracers()

plot_testcase()

advection_map()

advection_equatorial()

advection_error_convergence()
