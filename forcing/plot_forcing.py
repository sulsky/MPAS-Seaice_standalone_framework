from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import argparse
import sys
import numpy as np
from math import pi
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

def plot_forcing(filenameMesh, filenameData, filenameOut, varname, iTime):

    fileMesh = Dataset(filenameMesh,"r")

    nCells = len(fileMesh.dimensions["nCells"])
    nVertices = len(fileMesh.dimensions["nVertices"])

    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]

    verticesOnCell = fileMesh.variables["verticesOnCell"][:]

    latCell = fileMesh.variables["latCell"][:]
    lonCell = fileMesh.variables["lonCell"][:]

    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]

    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]

    fileMesh.close()

    verticesOnCell[:] -= 1

    for iCell in range(0,nCells):
        if (lonCell[iCell] < 0.0):
            lonCell[iCell] += 2.0*pi
    for iVertex in range(0,nVertices):
        if (lonVertex[iVertex] < 0.0):
            lonVertex[iVertex] += 2.0*pi

    xmin =  sys.float_info.max
    xmax = -sys.float_info.max
    ymin =  sys.float_info.max
    ymax = -sys.float_info.max
    for iVertex in range(0,nVertices):
        xmin = min(xmin,xVertex[iVertex])
        xmax = max(xmax,xVertex[iVertex])
        ymin = min(ymin,yVertex[iVertex])
        ymax = max(ymax,yVertex[iVertex])

    latmin =  sys.float_info.max
    latmax = -sys.float_info.max
    lonmin =  sys.float_info.max
    lonmax = -sys.float_info.max
    for iVertex in range(0,nVertices):
        latmin = min(latmin,latVertex[iVertex])
        latmax = max(latmax,latVertex[iVertex])
        lonmin = min(lonmin,lonVertex[iVertex])
        lonmax = max(lonmax,lonVertex[iVertex])
    fileData = Dataset(filenameData,"r")

    array = fileData.variables[varname][:]

    fileData.close()

    array = array[iTime,:]

    patchesNorth = []
    valuesNorth = []

    patchesSouth = []
    valuesSouth = []

    patchesLatLon = []
    valuesLatLon = []


    for iCell in range(0,nCells):

        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            if (lonVertex[iVertex] - lonCell[iCell] > 0.5*pi):
                vertices.append((lonVertex[iVertex]-2.0*pi,latVertex[iVertex]))
            elif (lonVertex[iVertex] - lonCell[iCell] < -0.5*pi):
                vertices.append((lonVertex[iVertex]+2.0*pi,latVertex[iVertex]))
            else:
                vertices.append((lonVertex[iVertex],latVertex[iVertex]))

        patchesLatLon.append(Polygon(vertices,closed=True))
        valuesLatLon.append(array[iCell])


        if (latCell[iCell] > 0.0):

            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append((xVertex[iVertex],yVertex[iVertex]))
            patchesNorth.append(Polygon(vertices,closed=True))
            valuesNorth.append(array[iCell])

        else:

            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append((yVertex[iVertex],xVertex[iVertex]))
            patchesSouth.append(Polygon(vertices,closed=True))
            valuesSouth.append(array[iCell])


    fig = plt.figure(figsize=(10, 7))
    gs = GridSpec(nrows=2, ncols=2)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[1, :])

    pcNorth = PatchCollection(patchesNorth, cmap="jet")
    pcNorth.set_array(np.array(valuesNorth))
    ax0.set_xlim((xmin, xmax))
    ax0.set_ylim((ymin, ymax))
    ax0.add_collection(pcNorth)
    ax0.set_aspect('equal')

    divider = make_axes_locatable(ax0)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pcNorth, cax=cax)

    pcSouth = PatchCollection(patchesSouth, cmap="jet")
    pcSouth.set_array(np.array(valuesSouth))
    ax1.set_xlim((xmin, xmax))
    ax1.set_ylim((ymin, ymax))
    ax1.add_collection(pcSouth)
    ax1.set_aspect('equal')

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pcSouth, cax=cax)

    pcLatLon = PatchCollection(patchesLatLon, cmap="jet")
    pcLatLon.set_array(np.array(valuesLatLon))
    ax2.set_xlim((lonmin, lonmax))
    ax2.set_ylim((latmin, latmax))
    ax2.add_collection(pcLatLon)
    ax2.set_aspect('equal')

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes('right', size='2.5%', pad=0.05)
    cb = fig.colorbar(pcLatLon, cax=cax)

    plt.tight_layout()
    plt.savefig(filenameOut)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m', dest='filenameMesh', help='')
    parser.add_argument('-i', dest='filenameData', help='')
    parser.add_argument('-o', dest='filenameOut', help='')
    parser.add_argument('-v', dest='varname', help='')
    parser.add_argument('-t', dest='iTime', type=int, default=0, help='')

    args = parser.parse_args()

    plot_forcing(args.filenameMesh, args.filenameData, args.filenameOut, args.varname, args.iTime)
