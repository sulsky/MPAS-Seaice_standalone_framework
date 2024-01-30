import os
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def check_results(n1,n2):

    try:
        import colorama
        colorama.init()
    except ImportError:
        pass

    print()
    message = "Check assembly results for %i and %i" %(n1,n2)
    print(message)
    print("="*len(message))

    cmd = "ncdiff -O output_%i/output.2000.nc output_%i/output.2000.nc diff.nc" %(n1,n2)
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

    message = "fieldAssemblyTest diff min/max: %g %g" %(fieldAssemblyTestMin, fieldAssemblyTestMax)
    try:
        import colorama
        print(colorama.Style.BRIGHT + colorama.Fore.MAGENTA + message + colorama.Style.RESET_ALL)
    except ImportError:
        print(message)
    print()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n1', dest="n1", type=int)
    parser.add_argument('--n2', dest="n2", type=int)
    args = parser.parse_args()

    check_results(args.n1, args.n2)
