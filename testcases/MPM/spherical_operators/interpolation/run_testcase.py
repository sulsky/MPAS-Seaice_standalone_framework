import sys
sys.path.append("../../../../utils/testcases/")
from get_testcase_data_spherical import get_testcase_data_spherical

from create_ic import create_ic
from create_particles import create_particles
from run_model import run_model
from interpolation_scaling import interpolation_scaling

#-------------------------------------------------------------------------------

def run_testcase():

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    create_particles()

    run_model()

    interpolation_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
