import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from create_ic import create_ic
from create_particles import create_particles
from run_model import run_model
from strain_map import strain_map
from strain_scaling import strain_scaling
from average_variational_stress import average_variational_stress
from stress_scaling import stress_scaling
from velocity_scaling import velocity_scaling
from velocity_map import velocity_map

#-------------------------------------------------------------------------------

def run_strain_testcase():

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    create_particles()

    run_model()

    average_variational_stress()

    strain_scaling()

    stress_scaling()

    velocity_scaling()

    velocity_map()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_strain_testcase()
