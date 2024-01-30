from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def empty_particle_file(filenameOut):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nParticles",0)

    fileOut.close()

#-------------------------------------------------------------------------------

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser()

    parser.add_argument('-o', dest="filenameOut", required=True)

    args = parser.parse_args()

    empty_particle_file(args.filenameOut)
