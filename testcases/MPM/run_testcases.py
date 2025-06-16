import sys

sys.path.append("../../utils/testcases")
from log_messages import log_message

import os
import subprocess
import time

#-------------------------------------------------------------------------------

def run_testcases():

    # List of test case directories
    testCaseDirs = [
        "advection",
        "advection_convergence",
        "assembly",
        "column_on_particles",
        "mpm_tracers",
        "polynya",
        "spherical_operators/interpolation",
        "spherical_operators/reconstruction",
        "standard_physics"
    ]
    testCaseDirs = [
        "column_on_particles"
    ]

    scriptDir = os.path.dirname(os.path.abspath(__file__))

    logFilename = "%s/log_mpm_testcases.txt" %(scriptDir)
    logFile = open(logFilename, "w")
    logFile.write("MPAS-Seaice-MPM test cases\n")
    logFile.write("##########################\n")

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        logFile.write("\nMPAS_SEAICE_EXECUTABLE not specified. Using standard location.\n")
    else:
        logFile.write("\nUsing MPAS-Seaice-MPM executable: %s\n" %(MPAS_SEAICE_EXECUTABLE))

    logFile.flush()
    logFile.close()

    for testDir in testCaseDirs:
        start_time = time.perf_counter()
        result = subprocess.run(["python", "run_testcase.py", "-l", logFilename], cwd=testDir)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        if result.returncode != 0:
            logFile = open(logFilename, "a")
            log_message("   Test case failed", "red", doPrint=False, logFile=logFile)
            logFile.flush()
            logFile.close()

        logFile = open(logFilename, "a")
        logFile.write(f"   Elapsed time: {elapsed_time:.4f} seconds\n")
        logFile.flush()
        logFile.close()

    logFile = open(logFilename, "a")
    logFile.write("\nTesting Completed\n")
    logFile.flush()
    logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcases()
