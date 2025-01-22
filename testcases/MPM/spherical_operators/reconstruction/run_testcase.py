from reconstruction_scaling import reconstruction_scaling
from reconstruction_map import reconstruction_map
import sys
sys.path.append("../interpolation")
from run_model import run_model
from create_ic import create_ic
from create_particles import create_particles
sys.path.append("../../../spherical_operators/strain_stress_divergence")
from get_testcase_data import get_testcase_data

#-------------------------------------------------------------------------------

def run_testcase():

    get_testcase_data()

    create_ic()

    create_particles()

    run_model()

    reconstruction_map()

    reconstruction_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
