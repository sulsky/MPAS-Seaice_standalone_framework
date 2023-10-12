from netCDF4 import Dataset
import math, os, sys
import random
import argparse

#--------------------------------------------------------------------

def create_ics(gridFilename):

    # grid in
    gridFile = Dataset(gridFilename, "r")

    nCells = len(gridFile.dimensions["nCells"])
    maxEdges = len(gridFile.dimensions["maxEdges"])

    gridFile.close()

    # ice out
    icFilename = "ic.nc"
    icFile = Dataset(icFilename, "w", format="NETCDF3_CLASSIC")

    icFile.createDimension("nCells", nCells)
    icFile.createDimension("maxEdges", maxEdges)

    fieldSubAssemblyTest = icFile.createVariable("fieldSubAssemblyTest", 'd', dimensions=("nCells","maxEdges"))

    for iCell in range(0,nCells):
        fieldSubAssemblyTest[iCell] = random.random()

    icFile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-m', dest="gridFilename")
    args = parser.parse_args()

    create_ics(args.gridFilename)
