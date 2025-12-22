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
import matplotlib
from matplotlib import colors
import argparse
from tqdm import tqdm

#-------------------------------------------------------------------------------

def iceberg_meltrates(filenameTemplate):

    # start plot
    fig, axis = plt.subplots()

    axis.set_facecolor('grey')


    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]

    fileMesh.close()

    filenames = sorted(glob.glob(filenameTemplate))

    icebergMeltRateCell = np.zeros(nCells)
    for filename in tqdm(filenames):
        fileIn = Dataset(filename,"r")

        icebergMeltRateCellIn = fileIn.variables["icebergMeltRateCell"][:,:]
        icebergMeltRateCell[:] += np.sum(icebergMeltRateCellIn, axis=0)

        fileIn.close()


    # plot mesh
    patches = []
    for iCell in range(0,nCells):
        if (latCell[iCell] < radians(-40.0)):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append([yVertex[iVertex],xVertex[iVertex]])
            patches.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))

    pc = PatchCollection(patches, match_original=True)

    axis.add_collection(pc)

    # plot melt rates
    patches = []
    colors = []
    for iCell in range(0,nCells):
        if (latCell[iCell] < radians(-40.0) and
            icebergMeltRateCell[iCell] > 0.0):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append([yVertex[iVertex],xVertex[iVertex]])
            patches.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))
            colors.append(icebergMeltRateCell[iCell])

    colors = np.array(colors)
    vmax = np.amax(colors)

    cmap = matplotlib.colormaps.get_cmap("jet")
    pc = PatchCollection(patches, match_original=True, cmap=cmap, norm=matplotlib.colors.LogNorm(vmin=vmax*0.001, vmax=vmax))
    pc.set_array(colors)

    axis.add_collection(pc)

    axis.autoscale_view()

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg melting")
    fig.colorbar(pc,label="Melt")

    plt.tight_layout()
    plt.savefig("iceberg_melting.png",dpi=600)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='')

    args = parser.parse_args()

    iceberg_meltrates(args.filenameTemplate)
