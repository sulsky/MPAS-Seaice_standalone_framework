from netCDF4 import Dataset, stringtochar
import geopandas as gpd
from math import degrees, sqrt, radians
import numpy as np
from shapely.geometry import Point
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import sys
import argparse
import os

#-------------------------------------------------------------------------------

def latlon_to_xyz(lat, lon):

    x = np.multiply(np.cos(lat),np.cos(lon))
    y = np.multiply(np.cos(lat),np.sin(lon))
    z = np.sin(lat)

    return x, y, z

#-------------------------------------------------------------------------------

def load_calving_data(calvingFilename):

    # Greene, C.A., Gardner, A.S., Schlegel, NJ. et al. Antarctic calving loss
    # rivals ice-shelf thinning. Nature 609, 948–953 (2022).
    # https://doi.org/10.1038/s41586-022-05037-w

    fileCalve = Dataset(calvingFilename,"r")

    nCalvingRates = len(fileCalve.dimensions["nRegions"])

    calvingNamesBytes = fileCalve.variables["regionName"][:]
    lonCalving = fileCalve.variables["lon"][:]
    latCalving = fileCalve.variables["lat"][:]
    calvingRate = fileCalve.variables["calvingRate"][:]

    fileCalve.close()

    calvingNames = []
    for iCalve in range(0,nCalvingRates):
        calvingNameStr = calvingNamesBytes[iCalve].tobytes().decode('utf-8')
        calvingNameStr = calvingNameStr.replace('\x00', '')
        calvingNames.append(calvingNameStr)
    calvingNames = np.array(calvingNames)

    # add some regions that help prevent transfer of calving mass across the peninsular
    # add West Graham Land
    calvingNames = np.append(calvingNames, "West Graham Land")
    lonCalving = np.append(lonCalving, -64.610219)
    latCalving = np.append(latCalving, -65.956736)
    calvingRate = np.append(calvingRate, 0.0)
    nCalvingRates += 1

    # add Eastern Graham Land
    calvingNames = np.append(calvingNames, "Eastern Graham Land")
    lonCalving = np.append(lonCalving, -58.378402)
    latCalving = np.append(latCalving, -63.898739)
    calvingRate = np.append(calvingRate, 0.0)
    nCalvingRates += 1

    return nCalvingRates, calvingNames, lonCalving, latCalving, calvingRate

#-------------------------------------------------------------------------------

def clean_ice_boundaries_region_names(nameIn, Subregion):

    # change region names to be consistent with calving names

    nameOut = nameIn.replace("_"," ")

    if (nameOut.find("Larsen") != -1):
        nameOut = nameOut.replace("Larsen","Larsen ")
        nameOut = nameOut.strip()
    elif (nameOut == "ClarkeBay"):
        nameOut = "Clarke Bay"
    elif (nameOut == "CapeWashington"):
        nameOut = "Cape Washington"
    elif (nameOut == "HarbordGlacier"):
        nameOut = "Harbord Glacier"
    elif (nameOut == "PourquoiPas"):
        nameOut = "Pourquoi Pas"
    elif (nameOut == "WattBay"):
        nameOut = "Watt Bay"
    elif (nameOut == "WilmaRobertDowner"):
        nameOut = "Wilma/Robert/Downer"
    elif (nameOut == "Fox"):
        if (Subregion == "East"):
            nameOut = "Fox Glacier"
        elif (Subregion == "West"):
            nameOut = "Fox Ice Stream"
        else:
            nameOut = "Not going to worry about this one"

    return nameOut

#-------------------------------------------------------------------------------

def ice_region_indices_corresponding_to_calving(nCalvingRates,
                                                calvingNames,
                                                iceRegions):

    # find the imbie regions corresponding to the calving regions

    regionIndices = []

    for iCalvingRates in range(0,nCalvingRates):

        calvingName = calvingNames[iCalvingRates]

        if (calvingName != "Other" and
            calvingName != "Antarctica"):

            nfound = 0
            for iRegion, iceRegionsRow in iceRegions.iterrows():
                regionName = iceRegionsRow['NAME']
                regionType = iceRegionsRow['TYPE']
                westEast = iceRegionsRow['Regions']

                regionName = clean_ice_boundaries_region_names(regionName, westEast)

                if (regionName == calvingName and
                    (regionType == "FL" or
                     (calvingName == "West Graham Land" and regionType == "GR") or
                     (calvingName == "Eastern Graham Land" and regionType == "GR"))):
                    regionIndices.append(iRegion)
                    nfound += 1

            if (nfound != 1):
                raise Exception("Could not find unique mapping from calving to region")

    return np.array(regionIndices)

#-------------------------------------------------------------------------------

def load_mpas_coastal_cells(meshFilename):

    # find the MPAS Antarctic coastal cells

    # load grid file data
    fileMesh = Dataset(meshFilename, "r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]
    lonCell = fileMesh.variables["lonCell"][:]

    areaCell = fileMesh.variables["areaCell"][:]

    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    cellsOnCell = fileMesh.variables["cellsOnCell"][:]-1

    fileMesh.close()

    coastalCellIndices = []
    latCellCoastal = []
    lonCellCoastal = []
    cellSizeCoastal = []

    for iCell in range(0,nCells):

        # boundary cell
        lBoundary = False
        for iCellOnCell in range(0,nEdgesOnCell[iCell]):
            iCell2 = cellsOnCell[iCell,iCellOnCell]
            if (iCell2 == -1):
                lBoundary = True

        if (lBoundary and
            degrees(latCell[iCell]) < -57.0):
            coastalCellIndices.append(iCell)
            latCellCoastal.append(latCell[iCell])
            lonCellCoastal.append(lonCell[iCell])
            cellSizeCoastal.append(sqrt(areaCell[iCell]))

    nCellsCoastal = len(coastalCellIndices)

    coastalCellIndices = np.array(coastalCellIndices)
    latCellCoastal = np.array(latCellCoastal)
    lonCellCoastal = np.array(lonCellCoastal)
    cellSizeCoastal = np.array(cellSizeCoastal)

    return nCellsCoastal, coastalCellIndices, latCellCoastal, lonCellCoastal, cellSizeCoastal

#-------------------------------------------------------------------------------

def plot_calving_cells(calvingCells,
                       iceRegions,
                       calvingPoints3031,
                       iceRegionsCalving,
                       coastalCellPoints3031):

    # plot
    ax = iceRegions.plot(figsize=(10, 10), edgecolor="grey", facecolor="lightblue", linewidth=0.2)

    # Add text labels at polygon centroids
    for idx, row in iceRegions.iterrows():
        centroid = row.geometry.centroid
        ax.text(
            centroid.x,
            centroid.y,
            str(row["NAME"]),  # replace with your actual column name
            fontsize=2,
            ha="center",
            va="center",
            color="black"
        )

    for iCalve in range(0,len(calvingCells)):
        plt.plot([calvingPoints3031.iloc[iCalve].x],
                 [calvingPoints3031.iloc[iCalve].y], marker="+", color="red", label="calving")
        plt.plot([iceRegionsCalving.iloc[iCalve].geometry.centroid.x],
                 [iceRegionsCalving.iloc[iCalve].geometry.centroid.y], marker="x", color="blue", label="regions")
        for iCell in list(calvingCells[iCalve]):
            plt.plot([coastalCellPoints3031.iloc[iCell].geometry.x],
                     [coastalCellPoints3031.iloc[iCell].geometry.y], marker="+", color="green", label="MPAS cell")
            plt.plot([calvingPoints3031.iloc[iCalve].x,coastalCellPoints3031.iloc[iCell].geometry.x],
                     [calvingPoints3031.iloc[iCalve].y,coastalCellPoints3031.iloc[iCell].geometry.y], color="magenta", lw=0.5)

    h,l = ax.get_legend_handles_labels()
    plt.legend(h[:3], l[:3])
    plt.title("MPAS cells, regions and calving linkage")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("calving_locations_antarctica.pdf")

#-------------------------------------------------------------------------------

def plot_calving_rate(nCells,
                      nEdgesOnCell,
                      verticesOnCell,
                      latCell,
                      xVertex,
                      yVertex,
                      calvingRateCells):

    # plot the calving rate
    patches = []
    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max
    colors = []
    for iCell in range(0,nCells):
        if (latCell[iCell] < radians(-50.0)):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                vertices.append([yVertex[iVertex],xVertex[iVertex]])
                xMin = min(xMin,xVertex[iVertex])
                xMax = max(xMax,xVertex[iVertex])
                yMin = min(yMin,yVertex[iVertex])
                yMax = max(yMax,yVertex[iVertex])
            patches.append(Polygon(vertices))
            colors.append(np.sum(calvingRateCells[iCell,:]))

    pc = PatchCollection(patches)
    pc.set_array(np.array(colors))

    fig, axis = plt.subplots(figsize=(10, 8))

    axis.add_collection(pc)

    axis.set_xlim((xMin,xMax))
    axis.set_ylim((yMin,yMax))

    axis.set_xlabel("x")
    axis.set_ylabel("y")

    axis.set_title("Antarctic calving rate")

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pc, cax=cax)

    cb.set_label("Calving rate (Gt/y)")

    plt.tight_layout()
    plt.savefig("mpas_calving_rate_antarctica.png",dpi=600)

#-------------------------------------------------------------------------------

def create_mpas_calving_file(calvingCells,
                             calvingRate,
                             calvingNames,
                             coastalCellIndices,
                             meshFilename,
                             calvingFilename,
                             multipleCalvingRegionsPerCell):

    nCalvingRegions = len(calvingCells)

    # create the output MPAS calving rate file

    fileMesh = Dataset(meshFilename, "r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]

    fileMesh.close()

    # number of calving regions per cell
    nCalvingRegionsPerCell = np.zeros(nCells,dtype="i")

    for iCalve in range(0,nCalvingRegions):
        nCellsCalve = len(calvingCells[iCalve])
        for iCell in list(calvingCells[iCalve]):
            nCalvingRegionsPerCell[coastalCellIndices[iCell]] += 1

    maxCalvingRegionsPerCell = np.amax(nCalvingRegionsPerCell)
    print("maxCalvingRegionsPerCell: ", maxCalvingRegionsPerCell)

    # calving rate per calving region per cell
    calvingRateCells = np.zeros((nCells,maxCalvingRegionsPerCell))
    calvingRegionsPerCell = np.zeros((nCells,maxCalvingRegionsPerCell),dtype="i")

    iCalvingRegionsPerCell = np.zeros(nCells,dtype="i")
    for iCalve in range(0,nCalvingRegions):
        nCellsCalve = len(calvingCells[iCalve])
        for iCell in list(calvingCells[iCalve]):
            calvingRateCells     [coastalCellIndices[iCell],iCalvingRegionsPerCell[coastalCellIndices[iCell]]] += calvingRate[iCalve] / float(nCellsCalve)
            calvingRegionsPerCell[coastalCellIndices[iCell],iCalvingRegionsPerCell[coastalCellIndices[iCell]]] = iCalve
            iCalvingRegionsPerCell[coastalCellIndices[iCell]] += 1

    # create file
    if (not multipleCalvingRegionsPerCell):
        maxCalvingRegionsPerCell = 1
        calvingRegionsPerCell = np.zeros((nCells,1))
        nCalvingRegions = 1
        calvingNames = ["NONE"]
        calvingRateCellsNew = np.zeros((nCells,1))
        calvingRateCellsNew[:,0] = np.sum(calvingRateCells,axis=1)
        calvingRateCells = calvingRateCellsNew
        calvingRate = np.array([np.sum(calvingRateCells)])
        nCalvingRegionsPerCell = np.clip(nCalvingRegionsPerCell,a_min=None,a_max=1)

    fileMPASCalving = Dataset(calvingFilename,"w",format="NETCDF3_CLASSIC")

    fileMPASCalving.totalCalvingRate = np.sum(calvingRateCells)
    fileMPASCalving.src = \
        "Greene, C.A., Gardner, A.S., Schlegel, NJ. et al. Antarctic calving loss " + \
        "rivals ice-shelf thinning. Nature 609, 948–953 (2022). " + \
        "https://doi.org/10.1038/s41586-022-05037-w"
    fileMPASCalving.regions = \
        "Mouginot, J., B. Scheuchl, and E. Rignot. 2017. MEaSUREs Antarctic Boundaries for IPY 2007-2009 " + \
        "from Satellite Radar, Version 2. [IceBoundaries_Antarctica_v02]. Boulder, Colorado USA. NASA National Snow " + \
        "and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/AXE4121732AD. [28th Oct 2025]"

    fileMPASCalving.createDimension("nCells", nCells)
    fileMPASCalving.createDimension("maxCalvingRegionsPerCell", maxCalvingRegionsPerCell)
    fileMPASCalving.createDimension("nCalvingRegions", nCalvingRegions)
    StrLen = 64
    fileMPASCalving.createDimension("StrLen", StrLen)

    var = fileMPASCalving.createVariable("nCalvingRegionsPerCell", "i", dimensions=["nCells"])
    var[:] = nCalvingRegionsPerCell[:]

    var = fileMPASCalving.createVariable("calvingRate", "d", dimensions=["nCells","maxCalvingRegionsPerCell"])
    var.units = "Gt/y"
    var[:] = calvingRateCells[:]

    var = fileMPASCalving.createVariable("calvingRegionIndex", "i", dimensions=["nCells","maxCalvingRegionsPerCell"])
    var[:] = calvingRegionsPerCell[:]

    var = fileMPASCalving.createVariable("calvingRegionNames", "c", dimensions=["nCalvingRegions","StrLen"])
    fixed = np.asarray(calvingNames, dtype=f'S{StrLen}')
    for iCalvingRegion in range(0,nCalvingRegions):
        var[:,:] = stringtochar(fixed)

    var = fileMPASCalving.createVariable("calvingRateRegions", "d", dimensions=["nCalvingRegions"])
    var.units = "Gt/y"
    var[:] = calvingRate[:]

    fileMPASCalving.close()

    plot_calving_rate(nCells,
                      nEdgesOnCell,
                      verticesOnCell,
                      latCell,
                      xVertex,
                      yVertex,
                      calvingRateCells)

#-------------------------------------------------------------------------------

def create_mpas_antarctic_calving_rates(meshFilename,
                                        calvingFilename,
                                        multipleCalvingRegionsPerCell):

    # get input calving and region files
    MPAS_SEAICE_STANDALONE_DATA = os.environ.get('MPAS_SEAICE_STANDALONE_DATA')
    if (MPAS_SEAICE_STANDALONE_DATA is None):
        raise Exception("MPAS_SEAICE_STANDALONE_DATA must be set")

    filenameCalvingRate = "%s/icebergs/Antarctica/calving_rates.nc" %(MPAS_SEAICE_STANDALONE_DATA)
    if (not os.path.isfile(filenameCalvingRate)):
        raise Exception("Could not find calving file: %s" %(filenameCalvingRate))

    filenameAntarcticRegions = "%s/icebergs/Antarctica/IceBoundaries_Antarctica_v02.shp" %(MPAS_SEAICE_STANDALONE_DATA)
    if (not os.path.isfile(filenameAntarcticRegions)):
        raise Exception("Could not find antarctic regions file: %s" %(filenameCalvingRate))

    # load MPAS coastal cells
    nCellsCoastal, coastalCellIndices, latCellCoastal, lonCellCoastal, cellSizeCoastal = \
        load_mpas_coastal_cells(meshFilename)
    print("nCellsCoastal: ", nCellsCoastal)

    # get calving rates
    nCalvingRates, calvingNames, lonCalving, latCalving, calvingRate = load_calving_data(filenameCalvingRate)
    print("nCalvingRates: ", nCalvingRates)

    # imbie region definitions
    # Mouginot, J., B. Scheuchl, and E. Rignot. 2017. MEaSUREs Antarctic Boundaries for IPY 2007-2009
    # from Satellite Radar, Version 2. [IceBoundaries_Antarctica_v02]. Boulder, Colorado USA. NASA National Snow
    # and Ice Data Center Distributed Active Archive Center. https://doi.org/10.5067/AXE4121732AD. [28th Oct 2025]
    iceRegions = gpd.read_file(filenameAntarcticRegions)
    print("nIceRegions: ", len(iceRegions))

    # find imbie region indices corresponding to the calving regions
    regionIndices = ice_region_indices_corresponding_to_calving(nCalvingRates,
                                                                calvingNames,
                                                                iceRegions)

    # subset of imbie regions corresponding to the calving regions
    iceRegionsCalving = iceRegions.iloc[regionIndices]
    iceRegionsCalving = iceRegionsCalving.reset_index(drop=True)
    print("nIceRegionsCalving: ", len(iceRegionsCalving))

    # Create a GeoDataFrame of MPAS coastal cells
    coastalCellCoords = []
    for iCell in range(0,nCellsCoastal):
        coastalCellCoords.append(Point(degrees(lonCellCoastal[iCell]),
                                       degrees(latCellCoastal[iCell])))

    coastalCellPoints = gpd.GeoDataFrame(
        geometry=coastalCellCoords,  # multiple points
        crs="EPSG:4326"
    )

    # Reproject both data frames to a common projected CRS - WGS 84 / Antarctic Polar Stereographic
    iceRegionsCalving3031 = iceRegionsCalving.to_crs("EPSG:3031")
    coastalCellPoints3031 = coastalCellPoints.to_crs("EPSG:3031")

    # calving locations as a geoseries
    calvingPointCoords = []

    for iCalve in range(0,regionIndices.shape[0]):
        calvingPointCoords.append(Point(lonCalving[iCalve], latCalving[iCalve]))

    calvingPoints = gpd.GeoSeries(
        calvingPointCoords,  # multiple points
        crs="EPSG:4326"
    )
    calvingPoints3031 = calvingPoints.to_crs("EPSG:3031")

    # get list of mpas cells associated with calving region
    calvingCells = []
    for iCalve in range(0,nCalvingRates):
        calvingCells.append(set())

    # nearest cell to calving region
    nearestCell = iceRegionsCalving3031.sjoin_nearest(coastalCellPoints3031,
                                                      how="left",
                                                      distance_col="distance")
    for iCalvingRegion, iceCalvingRegionsRow in nearestCell.iterrows():
        iCell = iceCalvingRegionsRow["index_right"]
        calvingCells[iCalvingRegion].add(iCell)

    # now all cells within a distance of calving region
    matches = []

    for iCell, coastalCellPointsRow in coastalCellPoints3031.iterrows():
        coastalCellGeometry = coastalCellPointsRow.geometry

        # Vectorized distance to all regions
        distances = iceRegionsCalving3031.distance(coastalCellGeometry)

        # Select only those within this point’s allowed distance
        nearbyRegions = iceRegionsCalving3031[distances < cellSizeCoastal[iCell]].copy()
        for iCalvingRegion in nearbyRegions.index.to_list():
            calvingCells[iCalvingRegion].add(iCell)

    plot_calving_cells(calvingCells,
                       iceRegions,
                       calvingPoints3031,
                       iceRegionsCalving,
                       coastalCellPoints3031)

    create_mpas_calving_file(calvingCells,
                             calvingRate,
                             calvingNames,
                             coastalCellIndices,
                             meshFilename,
                             calvingFilename,
                             multipleCalvingRegionsPerCell)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create MPAS-Seaice input iceberg calving file for Antarctica')

    parser.add_argument('-m', dest='meshFilename', required=True, help='Create MPAS-Seaice input iceberg calving file for Antarctica')
    parser.add_argument('-o', dest='calvingFilename', default="calving_mpas_antarctica.nc", help='MPAS calving input file name')
    parser.add_argument('-r', dest='multipleCalvingRegionsPerCell', action='store_true', help='List calving by source calving region')

    args = parser.parse_args()

    create_mpas_antarctic_calving_rates(args.meshFilename,
                                        args.calvingFilename,
                                        args.multipleCalvingRegionsPerCell)
