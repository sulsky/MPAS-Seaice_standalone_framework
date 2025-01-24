import os
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def check_results(n1, n2, method):

    cmd = "ncdiff -O output_%s_%i/output.2000.nc output_%s_%i/output.2000.nc diff.nc" %(method,n1,n2)
    os.system(cmd)

    cmd = "ncwa -O -y min diff.nc min.nc"
    os.system(cmd)

    cmd = "ncwa -O -y max diff.nc max.nc"
    os.system(cmd)

    filein = Dataset("min.nc","r")
    standardPhysicsTestMin = filein.variables["standardPhysicsTest"][:]
    filein.close()

    filein = Dataset("max.nc","r")
    standardPhysicsTestMax = filein.variables["standardPhysicsTest"][:]
    filein.close()

    print("standardPhysicsTest diff min/max: ", standardPhysicsTestMin, standardPhysicsTestMax)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n1', dest="n1", type=int)
    parser.add_argument('--n2', dest="n2", type=int)
    parser.add_argument('--o',  dest="operatorMethod", type=str)
    args = parser.parse_args()

    check_results(args.n1, args.n2, args.operatorMethod)
