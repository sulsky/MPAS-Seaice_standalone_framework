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

def iceberg_trajectories(filenameTemplate,
                         nskip):

    plt.rcParams["font.family"] = "Times New Roman"

    filenames = sorted(glob.glob(filenameTemplate))

    positions = {}

    vmin =  sys.float_info.max
    vmax = -sys.float_info.max

    for filename in tqdm(filenames):

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        if (nIcebergs > 0):

            icebergID = filein.variables["icebergID"][0,:]
            statusIB = filein.variables["statusIB"][0,:]
            posnIB = filein.variables["posnIBGeo"][0,:,:]
            icebergVolume = filein.variables["icebergVolume"][0,:]

            vmin = min(vmin,np.amin(icebergVolume))
            vmax = max(vmax,np.amax(icebergVolume))

            for iIceberg in range(0,nIcebergs):
                if (statusIB[iIceberg] == 1):
                    if (icebergID[iIceberg] not in positions):
                        positions[icebergID[iIceberg]] = {"x": [], "y": [], "v": []}
                    positions[icebergID[iIceberg]]["x"].append(posnIB[iIceberg,0])
                    positions[icebergID[iIceberg]]["y"].append(posnIB[iIceberg,1])
                    positions[icebergID[iIceberg]]["v"].append(icebergVolume[iIceberg])

        filein.close()


    # start plot
    fig, axis = plt.subplots()

    axis.set_facecolor('grey')


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
            lineSegments.append([[yVertex[iVertex1],xVertex[iVertex1]],
                                 [yVertex[iVertex2],xVertex[iVertex2]]])

    lc = LineCollection(lineSegments, color="black", linestyle='solid', linewidth=0.2)


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
    axis.add_collection(lc)

    # plot positions
    print("nIcebergs: ", len(positions))
    iIceberg = 0
    for icebergID, trajectory in tqdm(positions.items()):

        if (iIceberg % nskip == 0):
            lc = colored_line(trajectory["y"],
                              trajectory["x"],
                              trajectory["v"],
                              axis,
                              linewidth=0.2,
                              cmap="jet",
                              norm=colors.LogNorm(vmin=vmax*0.001, vmax=vmax))
        iIceberg += 1

    #axis.autoscale_view()
    axis.set_xlim(-3.5e6,3.5e6)
    axis.set_ylim(-3.5e6,3.5e6)

    axis.set_aspect("equal")
    axis.set_xlabel("x (m)")
    axis.set_ylabel("y (m)")
    axis.set_title("Iceberg trajectories")
    fig.colorbar(lc,label="Volume (m^3)")

    plt.tight_layout()
    plt.savefig("iceberg_trajectories.png",dpi=1200)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='')
    parser.add_argument('-n', dest='nskip', type=int, default=1, help='')

    args = parser.parse_args()

    iceberg_trajectories(args.filenameTemplate,
                         args.nskip)
