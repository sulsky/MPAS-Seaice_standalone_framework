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

def iceberg_meltrates():

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

    fileIn = Dataset("output/output.2000.nc")

    icebergMeltRateCell = fileIn.variables["icebergMeltRateCell"][-1,:]

    fileIn.close()


    # plot mesh
    patches = []
    colors = []
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
            colors.append(icebergMeltRateCell[iCell])

    pc = PatchCollection(patches, match_original=True, cmap="gist_stern")
    pc.set_array(colors)

    axis.add_collection(pc)

    axis.autoscale_view()

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg melting")
    fig.colorbar(pc,label="Melt")

    plt.tight_layout()
    plt.savefig("melting.png",dpi=600)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    iceberg_meltrates()
