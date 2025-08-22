from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def create_ics(gridFilename):

    uAirVelocity = 10.0
    vAirVelocity =  0.0

    # grid data
    fileGrid = Dataset(gridFilename,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    fileGrid.close()

    nCategories = 5

    iceAreaCategory = np.zeros((nCells,nCategories))
    iceVolumeCategory = np.zeros((nCells,nCategories))

    iceAreaCell = np.zeros((nCells))
    iceVolumeCell = np.zeros((nCells))

    for iCell in range(0, nCells):

        x = xCell[iCell]
        y = yCell[iCell]

        iceConcentration = 1
        iceThickness = 0.3

        iceAreaCategory[iCell,:] = 0.0
        iceAreaCategory[iCell,0] = iceConcentration
        iceAreaCell[iCell] = iceConcentration

        iceVolumeCategory[iCell,:] = 0.0
        iceVolumeCategory[iCell,0] = iceThickness
        iceVolumeCell[iCell] = iceThickness


    fileOut = Dataset("ic.nc", "w", format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("nCategories", nCategories)
    fileOut.createDimension("ONE", 1)

    var = fileOut.createVariable("iceAreaCategory", "d", dimensions=["nCells","nCategories","ONE"])
    var[:,:,0] = iceAreaCategory[:,:]

    var = fileOut.createVariable("iceVolumeCategory", "d", dimensions=["nCells","nCategories","ONE"])
    var[:,:,0] = iceVolumeCategory[:,:]

    var = fileOut.createVariable("iceAreaCell", "d", dimensions=["nCells"])
    var[:] = iceAreaCell[:]

    var = fileOut.createVariable("iceVolumeCell", "d", dimensions=["nCells"])
    var[:] = iceVolumeCell[:]

    fVertex = 0.0
    var = fileOut.createVariable("fVertex", "d", dimensions=["nVertices"])
    var[:] = fVertex

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
