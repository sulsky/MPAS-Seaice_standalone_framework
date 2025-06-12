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

    scriptDir = os.path.dirname(os.path.abspath(__file__))

    logFilename = "%s/log_mpm_testcases.txt" %(scriptDir)
    logFile = open(logFilename, "w")
    logFile.write("MPAS-Seaice-MPM test cases\n")
    logFile.write("##########################\n")
    logFile.flush()
    logFile.close()

    for testDir in testCaseDirs:
        start_time = time.perf_counter()
        result = subprocess.run(["python3", "run_testcase.py", "-l", logFilename], cwd=testDir)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time

        if result.returncode != 0:
            print(f"Test in {testDir} failed with return code {result.returncode}")

        logFile = open(logFilename, "a")
        logFile.write(f"   Elapsed time: {elapsed_time:.4f} seconds\n")
        logFile.flush()
        logFile.close()

    logFile = open(logFilename, "a")
    logFile.write("Testing Completed\n")
    logFile.flush()
    logFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcases()
