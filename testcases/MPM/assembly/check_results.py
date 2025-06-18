import sys

sys.path.append("../../../utils/testcases")
from log_messages import log_message

import os
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def check_results(n1,
                  n2,
                  logFile=None):

    try:
        import colorama
        colorama.init()
    except ImportError:
        pass

    print()
    message = "Check assembly results for %i and %i" %(n1,n2)
    print(message)
    print("="*len(message))

    filename1 = "output_%i/output.2000.nc" %(n1)
    filename2 = "output_%i/output.2000.nc" %(n2)
    if (not os.path.exists(filename1) or
        not os.path.exists(filename1)):
        log_message("   Missing check files", "red", logFile=logFile)
        raise Exception("Missing check files")
        return

    cmd = "ncdiff -O %s %s diff.nc" %(filename1, filename2)
    os.system(cmd)

    cmd = "ncwa -O -y min diff.nc min.nc"
    os.system(cmd)

    cmd = "ncwa -O -y max diff.nc max.nc"
    os.system(cmd)

    filein = Dataset("min.nc","r")
    fieldAssemblyTestMin = filein.variables["fieldAssemblyTest"][:]
    filein.close()

    filein = Dataset("max.nc","r")
    fieldAssemblyTestMax = filein.variables["fieldAssemblyTest"][:]
    filein.close()

    if (fieldAssemblyTestMin == 0 and
        fieldAssemblyTestMax == 0):
        message = "fieldAssemblyTest diff min/max: %g %g" %(fieldAssemblyTestMin, fieldAssemblyTestMax)
        log_message(message, "green", logFile=logFile)
    else:
        message = "fieldAssemblyTest diff min/max: %g %g" %(fieldAssemblyTestMin, fieldAssemblyTestMax)
        log_message(message, "red", logFile=logFile)

    print()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n1', dest="n1", type=int)
    parser.add_argument('--n2', dest="n2", type=int)
    args = parser.parse_args()

    check_results(args.n1, args.n2)
