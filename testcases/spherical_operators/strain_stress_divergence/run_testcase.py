import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from create_ic import create_ic
from run_model import run_model
from strain_stress_divergence_map import strain_stress_divergence_map
from strain_stress_divergence_scaling import strain_stress_divergence_scaling

#-------------------------------------------------------------------------------

def run_strain_stress_divergence_testcase():

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    run_model()

    strain_stress_divergence_map()

    strain_stress_divergence_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_strain_stress_divergence_testcase()
