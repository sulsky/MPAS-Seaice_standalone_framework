from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import numpy as np
import sys
from iceberg_plot_utils import projection, plot_limits, setup_maps_projection, mesh_patches
import argparse
from math import degrees

#-------------------------------------------------------------------------------

def initial_iceberg_locations(filenameTemplate,
                              location):

    filenames = sorted(glob.glob(filenameTemplate))

    src_crs, proj = setup_maps_projection(location)

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
    })

    initialPositions = {}

    vmin =  sys.float_info.max
    vmax = -sys.float_info.max

    for filename in filenames:

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        if (nIcebergs > 0):

            icebergID = filein.variables["icebergID"][0,:]
            statusIB = filein.variables["statusIB"][0,:]
            icebergVolume = filein.variables["icebergVolume"][0,:]
            latIceberg = filein.variables["latIceberg"][0,:]
            lonIceberg = filein.variables["lonIceberg"][0,:]

            vmin = min(vmin,np.amin(icebergVolume))
            vmax = max(vmax,np.amax(icebergVolume))

            for iIceberg in range(0,nIcebergs):
                if (statusIB[iIceberg] == 1):
                    if (icebergID[iIceberg] not in initialPositions):
                        initialPositions[icebergID[iIceberg]] = {"lat":latIceberg[iIceberg],
                                                                 "lon":lonIceberg[iIceberg]}


    filein.close()



    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])
    nEdges = len(fileMesh.dimensions["nEdges"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    cellsOnEdge = fileMesh.variables["cellsOnEdge"][:]-1
    verticesOnEdge = fileMesh.variables["verticesOnEdge"][:]-1
    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]
    latEdge = fileMesh.variables["latEdge"][:]

    fileMesh.close()


    # plot mesh
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

    # start plot
    fig = plt.figure(figsize=(5,6))
    axis = plt.axes(projection=proj)

    axis.set_facecolor('lightgrey')

    axis.add_collection(pc)
    axis.add_collection(lc)



    # plot initial locations
    xs = []
    ys = []
    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max
    for icebergID, position in initialPositions.items():
        x, y = projection(degrees(position["lat"]),
                          degrees(position["lon"]),
                          proj,src_crs)
        xs.append(x)
        ys.append(y)
        xMin = min(xMin,x)
        xMax = max(xMax,x)
        yMin = min(yMin,y)
        yMax = max(yMax,y)

    xMin, xMax, yMin, yMax = plot_limits(xMin, xMax, yMin, yMax)

    gl = axis.gridlines(linewidth=0.5,linestyle="dashed",draw_labels=False)

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

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='Input iceberg filename template to plot')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    initial_iceberg_locations(args.filenameTemplate,
                              args.location)
