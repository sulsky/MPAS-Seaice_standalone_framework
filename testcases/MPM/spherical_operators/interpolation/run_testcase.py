from create_ic import create_ic
from create_particles import create_particles
from run_model import run_model
from interpolation_scaling import interpolation_scaling
import sys
sys.path.append("../../../spherical_operators/strain_stress_divergence")
from get_testcase_data import get_testcase_data

#-------------------------------------------------------------------------------

def run_testcase():

    get_testcase_data()

    create_ic()

    create_particles()

    run_model()

    interpolation_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
