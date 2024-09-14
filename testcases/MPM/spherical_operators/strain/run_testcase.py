from get_testcase_data import get_testcase_data
from create_ic import create_ic
from create_particles import create_particles
from run_model import run_model
from strain_map import strain_map
from strain_scaling import strain_scaling
from stress_scaling import stress_scaling
from velocity_scaling import velocity_scaling

#-------------------------------------------------------------------------------

def run_strain_testcase():

    get_testcase_data()

    create_ic()

    create_particles()

    run_model()

    strain_scaling()

    stress_scaling()

    velocity_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_strain_testcase()
