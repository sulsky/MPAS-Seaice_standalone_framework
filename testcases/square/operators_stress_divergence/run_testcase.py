import argparse
from create_grids import create_grids
from create_ics import create_ics
from run_model import run_model
from stress_divergence_scaling import stress_divergence_scaling
from error_analysis_stress_divergence import error_analysis_stress_divergence

#-------------------------------------------------------------------------------

def run_testcase(meshType, meshScale, ignoreWeak):

    testName = create_grids(meshType, meshScale)

    create_ics(testName)

    run_model(testName, ignoreWeak)

    stress_divergence_scaling(testName, ignoreWeak)

    error_analysis_stress_divergence(testName, ignoreWeak)

#-------------------------------------------------------------------------------

if __name__ == "__main__" :

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='meshType', default="regular", choices=['regular','random','scale_x','scale_y'], help='')
    parser.add_argument('-s', dest='meshScale', type=float, help='')
    parser.add_argument('-w', dest='ignoreWeak', action='store_true', help='')

    args = parser.parse_args()

    run_testcase(args.meshType, args.meshScale, args.ignoreWeak)
