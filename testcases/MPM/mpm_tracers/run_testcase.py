import sys

sys.path.append("../../../testing")
from compare_mpas_files import compare_files
from testing_utils import get_domain, print_colour, create_test_directory, test_summary

sys.path.append("../../../utils/MPM/particle_initialization/")
from empty_particle_file import empty_particle_file

sys.path.append("../../../utils/testcases")
from log_messages import regression_summary

from run_model import run_model

import os
import argparse

#-------------------------------------------------------------------------------

def run_testcase(logFilenameOverview=None):

    if (logFilenameOverview is not None):
        logFileOverview = open(logFilenameOverview,"a")
        logFileOverview.write("\nMPM tracers test case\n")
        logFileOverview.write(  "=====================\n")
        logFileOverview.flush()

    # domains directory
    domainsDir = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
    if (domainsDir == None):
        raise Exception("Environment variable MPAS_SEAICE_DOMAINS_DIR must be set if no domains directory specified")
    if (not os.path.exists(domainsDir)):
        raise Exception("Requested domains directory does not exist")

    # get domain
    domain = "domain_QU120km"
    get_domain(domainsDir, domain)

    # empty particles
    empty_particle_file("particles.nc")

    # run models
    run_model()

    # check output
    # make a test directory
    testDir = "testDir"
    create_test_directory(testDir)
    os.chdir("../")

    # make log file
    logfile = open("log_test.txt", "w")
    title = "Test: mpm_tracers vs. mpas"
    print_colour(title, "title")
    logfile.write(title)

    # run comparison
    file1 = "output_mpm_tracers/output.2000.nc"
    file2 = "output_mpas/output.2000.nc"
    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile)
    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "mpm_tracers")
    regression_summary(nErrorsNonArray, nErrorsArray, logFileOverview, "mpm_tracers")

    if (os.path.isfile("vars_differ.nc")):
        cmd = "mv vars_differ.nc %s" %(testDir)
        os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='logFilename')
    args = parser.parse_args()

    run_testcase(args.logFilename)
