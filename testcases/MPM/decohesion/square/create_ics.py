from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def create_ics(gridFilename):

    uAir = 20.0
    vAir =  0.0

    # grid data
    fileGrid = Dataset(gridFilename,"r")

    nCells = len(fileGrid.dimensions["nCells"])

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    fileGrid.close()

    uAirVelocity = np.zeros(nCells)
    vAirVelocity = np.zeros(nCells)

    for iCell in range(0, nCells):

        x = xCell[iCell]
        y = yCell[iCell]
        uAirVelocity[iCell] = uAir
        vAirVelocity[iCell] = vAir

        if (x < 120000):
           uAirVelocity[iCell] = -uAir

    fileOut = Dataset("ic.nc", "w", format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells", nCells)

    var = fileOut.createVariable("uAirVelocity", "d", dimensions=["nCells"])
    var[:] = uAirVelocity

    var = fileOut.createVariable("vAirVelocity", "d", dimensions=["nCells"])
    var[:] = vAirVelocity

    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='gridFilename', required=True)

    args = parser.parse_args()

    create_ics(args.gridFilename)
