from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def fix_bathymetry(bathyFilename,
                   meshFilename,
                   correctionType):

    # mesh
    filemesh = Dataset("grid.nc","r")

    nCells = len(filemesh.dimensions["nCells"])

    nEdgesOnCell = filemesh.variables["nEdgesOnCell"][:]
    cellsOnCell  = filemesh.variables["cellsOnCell"][:]-1

    filemesh.close()

    # bathy file
    filein = Dataset("bathymetry.nc","r+")

    bathymetryVar = filein.variables["bathymetry"]

    bathymetry = bathymetryVar[:]

    nAvg              = np.zeros(nCells)
    bathymetryAvg     = np.zeros(nCells)
    bathymetryDeep    = np.zeros(nCells)
    bathymetryShallow = np.zeros(nCells)

    for iCell in range(0,nCells):

        if (bathymetry[iCell] > 0.0):

            bathymetryDeep   [iCell] = 0.0
            bathymetryShallow[iCell] = -1000.0

            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCellNeigh = cellsOnCell[iCell,iCellOnCell]

                print("  ", iCellOnCell, nEdgesOnCell[iCell], iCellNeigh, bathymetry[iCellNeigh])

                if (iCellNeigh >= 0 and
                    bathymetry[iCellNeigh] < 0.0):

                    # avg
                    bathymetryAvg[iCell] += bathymetry[iCellNeigh]
                    nAvg[iCell] += 1.0

                    # deep
                    bathymetryDeep   [iCell] = min(bathymetryDeep   [iCell],bathymetry[iCellNeigh])

                    # shallow
                    bathymetryShallow[iCell] = max(bathymetryShallow[iCell],bathymetry[iCellNeigh])

            if (nAvg[iCell] > 0):
                bathymetryAvg[iCell] /= nAvg[iCell]

            print(iCell,
                  bathymetry[iCell],
                  bathymetryAvg[iCell],
                  bathymetryDeep[iCell],
                  bathymetryShallow[iCell])

        else:

            bathymetryAvg    [iCell] = bathymetry[iCell]
            bathymetryDeep   [iCell] = bathymetry[iCell]
            bathymetryShallow[iCell] = bathymetry[iCell]

    if (correctionType == "avg"):
        bathymetryVar[:] = bathymetryAvg[:]
    elif (correctionType == "deep"):
        bathymetryVar[:] = bathymetryDeep[:]
    elif (correctionType == "shallow"):
        bathymetryVar[:] = bathymetryShallow[:]

    filein.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-b', dest='bathyFilename', required=True, help='')
    parser.add_argument('-m', dest='meshFilename',  required=True, help='')
    parser.add_argument('-t', dest='correctionType', choices=["avg","deep","shallow"], default="deep", help='')

    args = parser.parse_args()

    fix_bathymetry(args.bathyFilename,
                   args.meshFilename,
                   args.correctionType)
