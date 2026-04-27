from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import re
import numpy as np
import os
import sys
from math import radians
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
from iceberg_plot_utils import projection_scalar, projection_list, plot_limits
import argparse

#-------------------------------------------------------------------------------

def initial_iceberg_locations(location):

    filenames = sorted(glob.glob("./output/icebergs_output.*"))

    initialPositions = {}

    vmin =  sys.float_info.max
    vmax = -sys.float_info.max

    for filename in filenames:

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        if (nIcebergs > 0):

            icebergID = filein.variables["icebergID"][0,:]
            statusIB = filein.variables["statusIB"][0,:]
            posnIBInit = filein.variables["posnIBInit"][0,:,:]
            icebergVolume = filein.variables["icebergVolume"][0,:]

            vmin = min(vmin,np.amin(icebergVolume))
            vmax = max(vmax,np.amax(icebergVolume))

            for iIceberg in range(0,nIcebergs):
                if (statusIB[iIceberg] == 1):
                    if (icebergID[iIceberg] not in initialPositions):
                        initialPositions[icebergID[iIceberg]] = {"x":posnIBInit[iIceberg,0],
                                                                 "y":posnIBInit[iIceberg,1],
                                                                 "z":posnIBInit[iIceberg,2]}


    filein.close()

    # start plot
    fig, axis = plt.subplots()

    axis.set_facecolor('lightgrey')


    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]
    zVertex = fileMesh.variables["zVertex"][:]

    fileMesh.close()


    # plot mesh
    patches = []
    for iCell in range(0,nCells):
        if ((location == "antarctica" and latCell[iCell] < radians(-40.0)) or
            (location == "greenland"  and latCell[iCell] > radians( 40.0))):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                x, y = projection_scalar(xVertex[iVertex],
                                         yVertex[iVertex],
                                         zVertex[iVertex],
                                         location)
                vertices.append([x,y])
            patches.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))

    pc = PatchCollection(patches, match_original=True)

    axis.add_collection(pc)

    # plot initial locations
    xs = []
    ys = []
    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max
    for icebergID, position in initialPositions.items():
        x, y = projection_scalar(position["x"],
                                 position["y"],
                                 position["z"],
                                 location)
        xs.append(x)
        ys.append(y)
        xMin = min(xMin,x)
        xMax = max(xMax,x)
        yMin = min(yMin,y)
        yMax = max(yMax,y)

    xMin, xMax, yMin, yMax = plot_limits(xMin, xMax, yMin, yMax, enlarge=1.1)

    sc = axis.scatter(xs, ys, s=0.3, edgecolors='none')

    #axis.autoscale_view()
    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg initial positions")

    plt.tight_layout()
    plt.savefig("initial_positions.png",dpi=600)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    initial_iceberg_locations(args.location)
