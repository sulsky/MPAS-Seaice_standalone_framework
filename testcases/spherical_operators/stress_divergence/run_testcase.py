from get_testcase_data import get_testcase_data
from randomize_mesh import randomize_mesh
from create_ic import create_ic
from run_model import run_model
from stress_divergence_map import stress_divergence_map
from stress_divergence_scaling import stress_divergence_scaling
import argparse

#-------------------------------------------------------------------------------

def run_stress_divergence_testcase(meshType, meshScale):

    get_testcase_data()

    testName = randomize_mesh(meshType, meshScale)

    create_ic(testName)

    run_model(testName)

    stress_divergence_map(testName)

    stress_divergence_scaling(testName)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='meshType', default="regular", choices=['regular','random'], help='')
    parser.add_argument('-s', dest='meshScale', type=float, help='')

    args = parser.parse_args()

    run_stress_divergence_testcase(args.meshType, args.meshScale)
