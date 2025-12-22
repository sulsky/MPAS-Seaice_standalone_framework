from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def empty_particle_file(filenameOut,
                        decomposedDimensionName="nParticles"):

    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension(decomposedDimensionName,0)

    fileOut.close()

#-------------------------------------------------------------------------------

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser()

    parser.add_argument('-o', dest="filenameOut", required=True)
    parser.add_argument('-d', dest="decomposedDimensionName", default="nParticles")

    args = parser.parse_args()

    empty_particle_file(args.filenameOut,
                        args.decomposedDimensionName)
