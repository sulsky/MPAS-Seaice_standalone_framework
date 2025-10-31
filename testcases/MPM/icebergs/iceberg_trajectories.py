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

#-------------------------------------------------------------------------------

def iceberg_trajectories():

    filenames = sorted(glob.glob("./output/icebergs_output.*"))

    positions = {}

    for filename in filenames:

        print(filename)
        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        icebergID = filein.variables["icebergID"][0,:]
        posnIB = filein.variables["posnIB"][0,:,:]

        filein.close()

        for iIceberg in range(0,nIcebergs):

            if (icebergID[iIceberg] not in positions):
                positions[icebergID[iIceberg]] = {"x": [], "y": []}
            positions[icebergID[iIceberg]]["x"].append(posnIB[iIceberg,0])
            positions[icebergID[iIceberg]]["y"].append(posnIB[iIceberg,1])



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
        axis.plot(trajectory["y"], trajectory["x"], color="teal", linewidth=0.5)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg Trajectories")

    plt.tight_layout()
    plt.savefig("trajectories.png",dpi=600)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    iceberg_trajectories()
