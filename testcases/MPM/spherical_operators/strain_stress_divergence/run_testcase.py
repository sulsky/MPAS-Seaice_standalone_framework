from create_particles import create_particles
from run_model import run_model
from strain_stress_divergence_map import strain_stress_divergence_map
from strain_stress_divergence_scaling import strain_stress_divergence_scaling
import sys
sys.path.append("../../../spherical_operators/strain_stress_divergence")
from get_testcase_data import get_testcase_data
from create_ic import create_ic

#-------------------------------------------------------------------------------

def run_strain_stress_divergence_testcase():

    get_testcase_data()

    create_ic()

    create_particles()

    run_model()

    strain_stress_divergence_map()

    strain_stress_divergence_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_strain_stress_divergence_testcase()
