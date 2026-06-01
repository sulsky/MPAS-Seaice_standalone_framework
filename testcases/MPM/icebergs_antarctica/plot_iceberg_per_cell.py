from netCDF4 import Dataset
import argparse
import numpy as np
from iceberg_plot_utils import plot_cell_field

#-------------------------------------------------------------------------------

def per_cell_plots(filenameIn,
                   location):

    # load number of icebergs per cell
    fileCell = Dataset(filenameIn, "r")

    nIcebergs = len(fileCell.dimensions["nIcebergs"])

    globalCellIndexIB = fileCell.variables["globalCellIndexIB"][0,:]
    statusIB = fileCell.variables["statusIB"][0,:]
    icebergVolumeCell = fileCell.variables["icebergVolumeCell"][0,:]
    icebergAreaCell = fileCell.variables["icebergAreaCell"][0,:]
    posnIB = fileCell.variables["posnIBGeo"][0,:,:]

    fileCell.close()

    globalCellIndexIB = globalCellIndexIB.filled(-1)


    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])
    nEdges = len(fileMesh.dimensions["nEdges"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    cellsOnEdge = fileMesh.variables["cellsOnEdge"][:]-1
    verticesOnEdge = fileMesh.variables["verticesOnEdge"][:]-1
    indexToCellID = fileMesh.variables["indexToCellID"][:]
    areaCell = fileMesh.variables["areaCell"][:]
    latEdge = fileMesh.variables["latEdge"][:]
    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]

    fileMesh.close()

    globalCellIndexMap = {}
    for iCell in range(0,nCells):
        globalCellIndexMap[indexToCellID[iCell]] = iCell

    nIcebergsCell = np.zeros(nCells)
    for iIceberg in range(0,nIcebergs):
        if (statusIB[iIceberg] == 1):
            iCell = globalCellIndexMap[globalCellIndexIB[iIceberg]]
            nIcebergsCell[iCell] += 1

    # number plot plot
    plot_cell_field(nIcebergsCell,
                    "Iceberg number per cell",
                    r'Number',
                    "iceberg_num_per_cell.png",
                    location,
                    nEdges,
                    nCells,
                    nEdgesOnCell,
                    cellsOnEdge,
                    latEdge,
                    latCell,
                    verticesOnEdge,
                    verticesOnCell,
                    latVertex,
                    lonVertex)

    # volume plot plot
    plot_cell_field(icebergVolumeCell,
                    "Iceberg volume per cell",
                    r'Volume ($\mathrm{m}^3$)',
                    "iceberg_volume_per_cell.png",
                    location,
                    nEdges,
                    nCells,
                    nEdgesOnCell,
                    cellsOnEdge,
                    latEdge,
                    latCell,
                    verticesOnEdge,
                    verticesOnCell,
                    latVertex,
                    lonVertex)

    # area ratio plot
    icebergAreaCellRatio = icebergAreaCell[:] / areaCell[:]

    plot_cell_field(icebergAreaCellRatio,
                    "Iceberg area fraction",
                    r'Ratio',
                    "iceberg_area_fraction.png",
                    location,
                    nEdges,
                    nCells,
                    nEdgesOnCell,
                    cellsOnEdge,
                    latEdge,
                    latCell,
                    verticesOnEdge,
                    verticesOnCell,
                    latVertex,
                    lonVertex)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', dest='filenameIn', required=True, help='')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    per_cell_plots(args.filenameIn,
                   args.location)
