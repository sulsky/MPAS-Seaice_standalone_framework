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
import argparse
from tqdm import tqdm
from matplotlib import colors
from iceberg_plot_utils import projection_scalar, projection_list
from copy import copy

#-------------------------------------------------------------------------------

def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)

#-------------------------------------------------------------------------------

def plot_trajectories(location,
                      pc,
                      lc,
                      positions,
                      nskip,
                      fieldName,
                      vmin,
                      vmax,
                      colorbarTitle,
                      filenameOut,
                      xMin,
                      xMax,
                      yMin,
                      yMax):

    # start plot
    fig, axis = plt.subplots()

    axis.set_facecolor('grey')

    pcCopy = copy(pc)
    lcCopy = copy(lc)
    axis.add_collection(pcCopy)
    axis.add_collection(lcCopy)

    # plot positions
    iIceberg = 0
    for icebergID, trajectory in tqdm(positions.items()):

        if (iIceberg % nskip == 0):
            x, y = projection_list(trajectory["x"],
                                   trajectory["y"],
                                   trajectory["z"],
                                   location)
            lcLines = colored_line(x,
                                   y,
                                   trajectory[fieldName],
                                   axis,
                                   linewidth=0.1,
                                   cmap="jet",
                                   norm=colors.LogNorm(vmin=vmax*0.001, vmax=vmax))
        iIceberg += 1

    #axis.autoscale_view()
    axis.set_xlim(xMin,xMax)
    axis.set_ylim(yMin,yMax)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg trajectories")
    fig.colorbar(lcLines,label=colorbarTitle)

    plt.tight_layout()
    plt.savefig(filenameOut,dpi=600)
    plt.close()

#-------------------------------------------------------------------------------

def iceberg_trajectories(filenameTemplate,
                         nskip,
                         location):

    if (location == "antarctica"):
        xMin = -4e6
        xMax =  4e6
        yMin = -4e6
        yMax =  4e6
    elif (location == "greenland"):
        xMin = -1.5e6
        xMax =  1.5e6
        yMin = -2e6
        yMax =  1e6

    plt.rcParams["font.family"] = "Times New Roman"

    filenames = sorted(glob.glob(filenameTemplate))

    positions = {}

    vmin =  sys.float_info.max
    vmax = -sys.float_info.max
    smin =  sys.float_info.max
    smax = -sys.float_info.max

    for filename in tqdm(filenames):

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        if (nIcebergs > 0):

            icebergID = filein.variables["icebergID"][0,:]
            statusIB = filein.variables["statusIB"][0,:]
            posnIB = filein.variables["posnIBGeo"][0,:,:]
            icebergVolume = filein.variables["icebergVolume"][0,:]
            uVelocityIcebergGeo = filein.variables["uVelocityIcebergGeo"][0,:]
            vVelocityIcebergGeo = filein.variables["vVelocityIcebergGeo"][0,:]
            icebergSpeed = np.sqrt(np.add(np.power(uVelocityIcebergGeo,2),
                                          np.power(vVelocityIcebergGeo,2)))

            vmin = min(vmin,np.amin(icebergVolume))
            vmax = max(vmax,np.amax(icebergVolume))
            smin = min(smin,np.amin(icebergSpeed))
            smax = max(smax,np.amax(icebergSpeed))

            for iIceberg in range(0,nIcebergs):
                if (statusIB[iIceberg] == 1):
                    if (icebergID[iIceberg] not in positions):
                        positions[icebergID[iIceberg]] = {"x": [], "y": [], "z": [], "v": [], "s": []}
                    positions[icebergID[iIceberg]]["x"].append(posnIB[iIceberg,0])
                    positions[icebergID[iIceberg]]["y"].append(posnIB[iIceberg,1])
                    positions[icebergID[iIceberg]]["z"].append(posnIB[iIceberg,2])
                    positions[icebergID[iIceberg]]["v"].append(icebergVolume[iIceberg])
                    positions[icebergID[iIceberg]]["s"].append(icebergSpeed[iIceberg])

        filein.close()

    print("nIcebergs: ", len(positions))


    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])
    nEdges = len(fileMesh.dimensions["nEdges"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    cellsOnEdge = fileMesh.variables["cellsOnEdge"][:]-1
    verticesOnEdge = fileMesh.variables["verticesOnEdge"][:]-1
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]
    zVertex = fileMesh.variables["zVertex"][:]

    fileMesh.close()

    boundaryEdge = np.zeros(nEdges,dtype="i")
    for iEdge in range(0,nEdges):
        if (cellsOnEdge[iEdge,0] == -1 or
            cellsOnEdge[iEdge,1] == -1):
            boundaryEdge[iEdge] = 1

    lineSegments = []
    for iEdge in range(0,nEdges):
        if (boundaryEdge[iEdge] == 1):
            iVertex1 = verticesOnEdge[iEdge,0]
            iVertex2 = verticesOnEdge[iEdge,1]
            x1, y1 = projection_scalar(xVertex[iVertex1],
                                       yVertex[iVertex1],
                                       zVertex[iVertex1],
                                       location)
            x2, y2 = projection_scalar(xVertex[iVertex2],
                                       yVertex[iVertex2],
                                       zVertex[iVertex2],
                                       location)
            lineSegments.append([[x1,y2],
                                 [x1,y2]])

    lc = LineCollection(lineSegments, color="black", linestyle='solid', linewidth=0.2)

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


    plot_trajectories(location,
                      pc,
                      lc,
                      positions,
                      nskip,
                      "v",
                      vmin,
                      vmax,
                      "Volume (m^3)",
                      "iceberg_trajectories_volume.png",
                      xMin,
                      xMax,
                      yMin,
                      yMax)

    plot_trajectories(location,
                      pc,
                      lc,
                      positions,
                      nskip,
                      "s",
                      smin,
                      smax,
                      "Speed (m/s)",
                      "iceberg_trajectories_speed.png",
                      xMin,
                      xMax,
                      yMin,
                      yMax)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='')
    parser.add_argument('-n', dest='nskip', type=int, default=1, help='')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    iceberg_trajectories(args.filenameTemplate,
                         args.nskip,
                         args.location)
