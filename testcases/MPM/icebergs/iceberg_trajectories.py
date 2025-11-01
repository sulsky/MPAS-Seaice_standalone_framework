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

#-------------------------------------------------------------------------------

def iceberg_trajectories():

    filenames = sorted(glob.glob("./output/icebergs_output.*"))

    positions = {}

    vmin =  sys.float_info.max
    vmax = -sys.float_info.max

    for filename in filenames:

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        icebergID = filein.variables["icebergID"][0,:]
        posnIB = filein.variables["posnIB"][0,:,:]
        icebergVolume = filein.variables["icebergVolume"][0,:]

        vmin = min(vmin,np.amin(icebergVolume))
        vmax = max(vmax,np.amax(icebergVolume))

        filein.close()

        for iIceberg in range(0,nIcebergs):
            if (icebergID[iIceberg] not in positions):
                positions[icebergID[iIceberg]] = {"x": [], "y": [], "v": []}
            positions[icebergID[iIceberg]]["x"].append(posnIB[iIceberg,0])
            positions[icebergID[iIceberg]]["y"].append(posnIB[iIceberg,1])
            positions[icebergID[iIceberg]]["v"].append(icebergVolume[iIceberg])


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

    fileMesh.close()


    # plot mesh
    patches = []
    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max
    for iCell in range(0,nCells):
        if (latCell[iCell] < radians(-40.0)):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append([yVertex[iVertex],xVertex[iVertex]])
                xMin = min(xMin,xVertex[iVertex])
                xMax = max(xMax,xVertex[iVertex])
                yMin = min(yMin,yVertex[iVertex])
                yMax = max(yMax,yVertex[iVertex])
            patches.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))

    pc = PatchCollection(patches, match_original=True)

    axis.add_collection(pc)

    # plot trajectories
    for icebergID, trajectory in positions.items():

        sc = axis.scatter(trajectory["y"], trajectory["x"], c=trajectory["v"],
                          s=0.5, vmin=vmin, vmax=vmax, cmap="jet")

    axis.autoscale_view()

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg Trajectories")
    fig.colorbar(sc,label="Volume (m^3)")

    plt.tight_layout()
    plt.savefig("trajectories.png",dpi=600)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    iceberg_trajectories()
