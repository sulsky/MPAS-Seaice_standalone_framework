import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_square import get_testcase_data_square

sys.path.append("../../../square/square_quadhex")
from create_grids import create_grids

from create_ics import create_ics
from create_particles import create_particles
from run_model import run_model
from set_difference_fields import set_difference_fields
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_square_quadhex_testcase():

    print("Create grids")
    print("=================")
    create_grids()

    print("Create ics")
    print("=================")
    create_ics()

    print("Create particles")
    print("=================")
    create_particles()

    print("Run model")
    print("=================")
    run_model()

    print("Set difference fields")
    print("=================")
    set_difference_fields()

    print("Plot testcase")
    print("=================")
    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_square_quadhex_testcase()
