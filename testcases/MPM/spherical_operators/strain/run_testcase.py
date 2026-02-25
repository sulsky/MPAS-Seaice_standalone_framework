import sys

sys.path.append("../../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from create_ic import create_ic
from create_particles import create_particles
from run_model import run_model
from strain_map import strain_map
from strain_scaling import strain_scaling
from average_variational_stress import average_variational_stress
from stress_scaling import stress_scaling
from stress_map import stress_map
from velocity_scaling import velocity_scaling
from velocity_map import velocity_map
import argparse

#-------------------------------------------------------------------------------

def run_strain_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nstrain test case\n")
        logFile.write(  "================\n")
        logFile.flush()
    else:
        logFile = None

    get_testcase_data_spherical(getGraphFiles=False)

    print("Create ICs...")
    print("=================")
    create_ic()

    print("Create particles...")
    print("=================")
    create_particles()

    print("Run model...")
    print("=================")
    run_model(logFile)

    print("Velocity scaling...")
    print("=================")
    velocity_scaling()

    print("Velocity map...")
    print("=================")
    velocity_map()

    print("Strain scaling...")
    print("=================")
    strain_scaling()

    print("Strain map...")
    print("=================")
    strain_map()

    print("Average variational stress...")
    print("=================")
    average_variational_stress()

    print("Stress scaling...")
    print("=================")
    stress_scaling()

    print("Stress map...")
    print("=================")
    stress_map()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='logFilename')
    args = parser.parse_args()

    run_strain_testcase(args.logFilename)
