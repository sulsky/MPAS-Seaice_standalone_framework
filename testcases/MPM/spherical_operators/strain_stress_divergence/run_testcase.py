import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

sys.path.append("../../../spherical_operators/strain_stress_divergence")
from create_ic import create_ic

from create_particles import create_particles
from run_model import run_model
from strain_stress_divergence_map import strain_stress_divergence_map
from strain_stress_divergence_scaling import strain_stress_divergence_scaling

import argparse

#-------------------------------------------------------------------------------

def run_strain_stress_divergence_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nstrain stress divergence test case\n")
        logFile.write(  "==================================\n")
        logFile.flush()
    else:
        logFile = None

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    create_particles()

    run_model(logFile)

    strain_stress_divergence_map()

    strain_stress_divergence_scaling()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='logFilename')
    args = parser.parse_args()

    run_strain_stress_divergence_testcase(args.logFilename)
