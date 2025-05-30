import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from create_ic import create_ic
from run_model import run_model
from average_variational_strains import average_variational_strains
from average_weak_strains import average_weak_strains
from strain_map import strain_map
from strain_scaling import strain_scaling

#-------------------------------------------------------------------------------

def run_strain_testcase():

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    run_model()

    average_variational_strains()

    average_weak_strains()

    strain_map()

    strain_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_strain_testcase()
