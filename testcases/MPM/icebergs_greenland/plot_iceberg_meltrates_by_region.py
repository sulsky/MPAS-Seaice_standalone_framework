from netCDF4 import Dataset, chartostring
import matplotlib.pyplot as plt
import glob
import numpy as np
import argparse
from tqdm import tqdm
from iceberg_plot_utils import plot_limits, setup_maps_projection, plot_cell_field_axis, data_range

#-------------------------------------------------------------------------------

def iceberg_meltrates_by_region(filenameTemplate,
                                location):

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
    areaCell = fileMesh.variables["areaCell"][:]
    latVertex = np.degrees(fileMesh.variables["latVertex"][:])
    lonVertex = np.degrees(fileMesh.variables["lonVertex"][:])
    latEdge = fileMesh.variables["latEdge"][:]

    fileMesh.close()

    filenames = sorted(glob.glob(filenameTemplate))

    fileIn = Dataset(filenames[0],"r")
    calvingRegionNames = chartostring(fileIn.variables["calvingRegionNames"][:])
    fileIn.close()



    icebergMeltRateCellByRegion = {}

    icebergMeltRateCell = np.zeros(nCells)
    nAvg = 0
    for filename in tqdm(filenames):
        fileIn = Dataset(filename,"r")

        icebergMeltRateCellIn = fileIn.variables["icebergMeltRateCell"][:,:]
        icebergMeltRateCell[:] += np.sum(icebergMeltRateCellIn, axis=0)

        fileIn.close()

        nAvg += icebergMeltRateCellIn.shape[0]

    icebergMeltRateCell[:] /= nAvg

    icebergMeltRateCellByRegion["All"] = {"data":icebergMeltRateCell,
                                          "name":"All"}


    nRegions = 7
    for iRegion in range(0,nRegions):

        icebergMeltRateCell = np.zeros(nCells)
        nAvg = 0
        for filename in tqdm(filenames):
            fileIn = Dataset(filename,"r")

            icebergMeltRateCellIn = fileIn.variables["icebergMeltRateCellByRegion"][:,:,iRegion]
            icebergMeltRateCell[:] += np.sum(icebergMeltRateCellIn, axis=0)

            fileIn.close()

            nAvg += icebergMeltRateCellIn.shape[0]

        icebergMeltRateCell[:] /= nAvg

        icebergMeltRateCellByRegion[calvingRegionNames[iRegion]] = {"data":icebergMeltRateCell,
                                                                    "name":calvingRegionNames[iRegion]}


    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
    })

    src_crs, proj = setup_maps_projection(location)

    xMin, xMax, yMin, yMax = data_range(icebergMeltRateCellByRegion["All"]["data"],
                                        src_crs,
                                        proj,
                                        location,
                                        nCells,
                                        nEdgesOnCell,
                                        verticesOnCell,
                                        latCell,
                                        latVertex,
                                        lonVertex)
    xMin, xMax, yMin, yMax = plot_limits(xMin, xMax, yMin, yMax)

    regions = ["CE","CW","NE","NO","NW","SE","SW","All"]

    fig, axes = plt.subplots(2,4, subplot_kw={'projection':proj}, figsize=(10,6))

    for iRegion in tqdm(range(0,len(regions))):
        i = iRegion // 4
        j = iRegion % 4

        plot_cell_field_axis(fig,
                             axes[i,j],
                             icebergMeltRateCellByRegion[regions[iRegion]]["data"],
                             regions[iRegion],
                             r'Melt flux ($\mathrm{kg}/\mathrm{m}^2/\mathrm{s}$)',
                             src_crs,
                             proj,
                             location,
                             nEdges,
                             nCells,
                             nEdgesOnCell,
                             cellsOnEdge,
                             latEdge,
                             latCell,
                             verticesOnEdge,
                             verticesOnCell,
                             latVertex,
                             lonVertex,
                             xMinIn=xMin,
                             xMaxIn=xMax,
                             yMinIn=yMin,
                             yMaxIn=yMax)

    plt.tight_layout()
    plt.savefig("iceberg_meltrate_by_region.png",dpi=300)
    plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', required=True, help='')
    parser.add_argument('-l', dest='location', choices=["antarctica","greenland"], default="antarctica", help='')

    args = parser.parse_args()

    iceberg_meltrates_by_region(args.filenameTemplate,
                                args.location)
