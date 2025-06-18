
from log_messages import log_message
import os
import subprocess
import time

#-------------------------------------------------------------------------------

def execute_model(MPAS_SEAICE_EXECUTABLE,
                  nProcs=1,
                  logFile=None):

    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND == "mpirun"):
        cmd = ["mpirun", "-oversubscribe", "-np", str(nProcs), MPAS_SEAICE_EXECUTABLE]
    elif (MPAS_SEAICE_TESTCASES_RUN_COMMAND == "srun"):
        cmd = ["srun", "--nodes=1", "--cpus-per-task=1", "--ntasks-per-node=%i" %(nProcs), MPAS_SEAICE_EXECUTABLE]
    else:
        raise Exception("Unsupported MPAS_SEAICE_TESTCASES_RUN_COMMAND type: %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND))

    print(' '.join(cmd))
    start = time.perf_counter()
    result = subprocess.run(cmd)
    if (logFile is not None and
        result.returncode != 0):
        log_message("   ERROR: Executable did not complete.", "red", logFile=logFile)
    end = time.perf_counter()
    print(f"Elapsed time: {end - start:.6f} seconds")

#-------------------------------------------------------------------------------
