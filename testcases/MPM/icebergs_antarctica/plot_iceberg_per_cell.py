from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse
import sys
import numpy as np
from iceberg_plot_utils import projection_scalar, projection_list, plot_limits

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

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]
    zVertex = fileMesh.variables["zVertex"][:]
    indexToCellID = fileMesh.variables["indexToCellID"][:]
    areaCell = fileMesh.variables["areaCell"][:]

    fileMesh.close()

    globalCellIndexMap = {}
    for iCell in range(0,nCells):
        globalCellIndexMap[indexToCellID[iCell]] = iCell

    nIcebergsCell = np.zeros(nCells)
    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max
    for iIceberg in range(0,nIcebergs):
        if (statusIB[iIceberg] == 1):
            iCell = globalCellIndexMap[globalCellIndexIB[iIceberg]]
            nIcebergsCell[iCell] += 1
            x, y = projection_scalar(posnIB[iIceberg,0],
                                     posnIB[iIceberg,1],
                                     posnIB[iIceberg,2],
                                     location)
            xMin = min(xMin,x)
            xMax = max(xMax,x)
            yMin = min(yMin,y)
            yMax = max(yMax,y)


    # plot mesh
    patchesCell = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            x, y = projection_scalar(xVertex[iVertex],
                                     yVertex[iVertex],
                                     zVertex[iVertex],
                                     location)

            vertices.append([x,y])
        patchesCell.append(Polygon(vertices, closed=True))

    xMin, xMax, yMin, yMax = plot_limits(xMin, xMax, yMin, yMax, enlarge=1.1)

    # number plot plot
    fig, axis = plt.subplots()

    axis.set_facecolor('grey')

    pcCell = PatchCollection(patchesCell, match_original=True, cmap=plt.get_cmap('jet'))
    pcCell.set_array(nIcebergsCell)
    axis.add_collection(pcCell)

    #axis.autoscale_view()
    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg number per cell")
    fig.colorbar(pcCell,label="num icebergs")

    plt.tight_layout()
    plt.savefig("iceberg_num_per_cell.png",dpi=1200)

    # volume plot plot
    fig, axis = plt.subplots()

    axis.set_facecolor('grey')

    pcCell = PatchCollection(patchesCell, match_original=True, cmap=plt.get_cmap('jet'))
    pcCell.set_array(icebergVolumeCell)
    axis.add_collection(pcCell)

    #axis.autoscale_view()
    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg volume per cell")
    fig.colorbar(pcCell,label="m^3")

    plt.tight_layout()
    plt.savefig("iceberg_volume_per_cell.png",dpi=1200)

    # area ratio plot
    icebergAreaCellRatio = icebergAreaCell[:] / areaCell[:]

    fig, axis = plt.subplots()

    axis.set_facecolor('grey')

    pcCell = PatchCollection(patchesCell, match_original=True, cmap=plt.get_cmap('jet'))
    pcCell.set_array(icebergAreaCellRatio)
    axis.add_collection(pcCell)

    #axis.autoscale_view()
    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg area fraction")
    fig.colorbar(pcCell,label="-")

    plt.tight_layout()
    plt.savefig("iceberg_area_fraction.png",dpi=1200)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', dest='filenameIn', required=True, help='')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    per_cell_plots(args.filenameIn,
                   args.location)
