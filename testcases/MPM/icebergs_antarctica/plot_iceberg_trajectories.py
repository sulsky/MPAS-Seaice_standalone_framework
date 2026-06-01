from netCDF4 import Dataset, chartostring
import matplotlib.pyplot as plt
import glob
import numpy as np
import sys
from matplotlib.collections import LineCollection
import argparse
from tqdm import tqdm
from matplotlib import colors
from iceberg_plot_utils import plot_limits, projection, projection_list, setup_maps_projection, mesh_patches
from copy import copy

secsToDays = 1.0 / 86400.0

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
                      proj,
                      src_crs,
                      pc,
                      lc,
                      positions,
                      nskip,
                      fieldName,
                      vmin,
                      vmax,
                      colorbarTitle,
                      filenameOut,
                      calvingRegionIndices=None,
                      useLogPlot=True):

    # start plot
    fig = plt.figure(figsize=(7,6))
    axis = plt.axes(projection=proj)

    axis.set_facecolor('grey')

    pcCopy = copy(pc)
    lcCopy = copy(lc)
    axis.add_collection(pcCopy)
    axis.add_collection(lcCopy)


    if (useLogPlot):
        norm = colors.LogNorm(vmin=vmax*0.001, vmax=vmax)
    else:
        norm = colors.Normalize(vmin=vmin, vmax=730)

    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max
    iIceberg = 0
    for icebergID, trajectory in tqdm(positions.items()):
        if (calvingRegionIndices is None or
            trajectory["o"][0] in calvingRegionIndices):

            if (iIceberg % nskip == 0):
                x, y = projection_list(np.degrees(trajectory["lat"]),
                                       np.degrees(trajectory["lon"]),
                                       proj,src_crs)
                lcLines = colored_line(x,
                                       y,
                                       trajectory[fieldName],
                                       axis,
                                       linewidth=0.1,
                                       cmap="viridis",
                                       norm=norm)
                xMin = min(xMin,np.amin(x))
                xMax = max(xMax,np.amax(x))
                yMin = min(yMin,np.amin(y))
                yMax = max(yMax,np.amax(y))
        iIceberg += 1

    xMin, xMax, yMin, yMax = plot_limits(xMin, xMax, yMin, yMax)

    gl = axis.gridlines(linewidth=0.5,linestyle="dashed",draw_labels=False)

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
                         location,
                         calvingRegion):

    src_crs, proj = setup_maps_projection(location)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
    })

    filenames = sorted(glob.glob(filenameTemplate))

    positions = {}

    vmin =  sys.float_info.max
    vmax = -sys.float_info.max
    smin =  sys.float_info.max
    smax = -sys.float_info.max
    amin =  sys.float_info.max
    amax = -sys.float_info.max

    if (calvingRegion is not None):
        filein = Dataset(filenames[0],"r")
        nCalvingRegions = len(filein.dimensions["nCalvingRegions"])
        calvingRegionNames = filein.variables["calvingRegionNames"][:]
        calvingRegionNames = chartostring(calvingRegionNames)
        filein.close()
        calvingRegionIndices = []
        for iName in range(0,nCalvingRegions):
            if (calvingRegionNames[iName] == calvingRegion):
                calvingRegionIndices.append(iName)
    else:
        calvingRegionIndices = None

    for filename in tqdm(filenames):

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        if (nIcebergs > 0):

            icebergID = filein.variables["icebergID"][0,:]
            statusIB = filein.variables["statusIB"][0,:]
            icebergCalvingRegionIndex = filein.variables["icebergCalvingRegionIndex"][0,:]
            latIceberg = filein.variables["latIceberg"][0,:]
            lonIceberg = filein.variables["lonIceberg"][0,:]
            icebergVolume = filein.variables["icebergVolume"][0,:]
            icebergSpeed = filein.variables["icebergSpeed"][0,:]
            icebergAge = filein.variables["icebergAge"][0,:]
            icebergAge[:] *= secsToDays

            vmin = min(vmin,np.amin(icebergVolume))
            vmax = max(vmax,np.amax(icebergVolume))
            smin = min(smin,np.amin(icebergSpeed))
            smax = max(smax,np.amax(icebergSpeed))
            amin = min(amin,np.amin(icebergAge))
            amax = max(amax,np.amax(icebergAge))

            for iIceberg in range(0,nIcebergs):
                if (statusIB[iIceberg] == 1):
                    if (icebergID[iIceberg] not in positions):
                        positions[icebergID[iIceberg]] = {"o": [],
                                                          "lat": [],
                                                          "lon": [],
                                                          "v": [],
                                                          "s": [],
                                                          "a": []}
                    positions[icebergID[iIceberg]]["o"].append(icebergCalvingRegionIndex[iIceberg])
                    positions[icebergID[iIceberg]]["lat"].append(latIceberg[iIceberg])
                    positions[icebergID[iIceberg]]["lon"].append(lonIceberg[iIceberg])
                    positions[icebergID[iIceberg]]["v"].append(icebergVolume[iIceberg])
                    positions[icebergID[iIceberg]]["s"].append(icebergSpeed[iIceberg])
                    positions[icebergID[iIceberg]]["a"].append(icebergAge[iIceberg])

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
    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]
    latEdge = fileMesh.variables["latEdge"][:]
    lonEdge = fileMesh.variables["lonEdge"][:]

    fileMesh.close()

    pc, lc = mesh_patches(location,
                          proj,
                          src_crs,
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

    plot_trajectories(location,
                      proj,
                      src_crs,
                      pc,
                      lc,
                      positions,
                      nskip,
                      "v",
                      vmin,
                      vmax,
                      r'Volume ($\mathrm{m}^3$)',
                      "iceberg_trajectories_volume.png",
                      calvingRegionIndices=calvingRegionIndices)

    plot_trajectories(location,
                      proj,
                      src_crs,
                      pc,
                      lc,
                      positions,
                      nskip,
                      "s",
                      smin,
                      smax,
                      r'Speed ($\mathrm{m}/\mathrm{s}$)',
                      "iceberg_trajectories_speed.png",
                      calvingRegionIndices=calvingRegionIndices)

    plot_trajectories(location,
                      proj,
                      src_crs,
                      pc,
                      lc,
                      positions,
                      nskip,
                      "a",
                      amin,
                      amax,
                      r'Age (days)',
                      "iceberg_trajectories_age.png",
                      calvingRegionIndices=calvingRegionIndices,
                      useLogPlot=False)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='Input iceberg filename template to plot')
    parser.add_argument('-n', dest='nskip', type=int, default=1, help='Stride for plotting iceberg subset')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='Plotting location')
    parser.add_argument('-c', dest='calvingRegion', default=None, help='plots only icebergs from this calving region')

    args = parser.parse_args()

    iceberg_trajectories(args.filenameTemplate,
                         args.nskip,
                         args.location,
                         args.calvingRegion)
