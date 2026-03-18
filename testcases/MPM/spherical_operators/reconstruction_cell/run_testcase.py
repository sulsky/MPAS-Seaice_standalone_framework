import sys

sys.path.append("../../../../utils/testcases/")
from get_testcase_data_spherical import get_testcase_data_spherical
from log_messages import log_message

from create_ic import create_ic
from create_particles import create_particles
from run_model import run_model
from reconstruction_cell_scaling import reconstruction_cell_scaling

import os
import argparse

#-------------------------------------------------------------------------------

def run_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nReconstruction cell test case\n")
        logFile.write(  "=======================\n")
        logFile.flush()
    else:
        logFile = None

    get_testcase_data_spherical(getGraphFiles=False)

    create_ic()

    create_particles()

    run_model(logFile)

    reconstruction_cell_scaling()

    if (logFile is not None):
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        log_message("Check plots", "yellow", logFile=logFile)
        log_message("  "+scriptDir+"/reconstruction_cell_scaling.png", "magenta", logFile=logFile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='logFilename')
    args = parser.parse_args()

    run_testcase(args.logFilename)
