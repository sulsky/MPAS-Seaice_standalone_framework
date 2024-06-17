from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

#---------------------------------------------------------

def get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, array):

    # find initial point
    radius = 6371229.0
    xStart = radius
    yStart = 0.0
    zStart = 0.0

    minDistance = 1e30
    iCellStart = -1
    for iCell in range(0,nCells):

        distance = math.sqrt(math.pow((xCell[iCell] - xStart),2) + \
                             math.pow((yCell[iCell] - yStart),2) + \
                             math.pow((zCell[iCell] - zStart),2))

        if (distance < minDistance):
            minDistance = distance
            iCellStart = iCell


    iCell = iCellStart

    lats = []
    lons = []
    y = []
    deltaLon = 0.0

    while True:

        lats.append(latCell[iCell])
        lons.append(lonCell[iCell])
        y.append(array[iCell])

        iCellNext = -1
        nextSmallestLatitude = 1e30

        for iCellOnCell in range(0,nEdgesOnCell[iCell]):

            iCellNeighbour = cellsOnCell[iCell,iCellOnCell] - 1

            lonNeighbour = lonCell[iCellNeighbour]

            if (lonNeighbour < lonCell[iCell] - math.pi):
                lonNeighbour = lonNeighbour + 2.0 * math.pi
            if (lonNeighbour > lonCell[iCell] + math.pi):
                lonNeighbour = lonNeighbour - 2.0 * math.pi

            if (lonNeighbour > lonCell[iCell]):

                if (math.fabs(latCell[iCellNeighbour]) < nextSmallestLatitude):

                    nextSmallestLatitude = math.fabs(latCell[iCellNeighbour])
                    iCellNext = iCellNeighbour
                    lonNeighbourAll = lonNeighbour

        deltaLon = deltaLon + (lonNeighbourAll - lonCell[iCell])

        if (deltaLon > 2.0 * math.pi):
            break

        iCell = iCellNext

    points = zip(lons,y)
    sorted_points = sorted(points)
    lons = [point[0] for point in sorted_points]
    y    = [point[1] for point in sorted_points]

    return lons, y

#---------------------------------------------------------

def advection_equatorial():

    #res = "2562"
    #res = "10242"
    res = "40962"
    #res = "163842"

    experiment1 = "cosine_bell"
    experiment2 = "slotted_cylinder"

    iTime = -1

    # cosine bell
    filename = "./output_%s_%s/output.2000.nc" %(experiment1,res)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    nCells = len(fileMPAS.dimensions["nCells"])

    latCell = fileMPAS.variables["latCell"][:]
    lonCell = fileMPAS.variables["lonCell"][:]

    xCell = fileMPAS.variables["xCell"][:]
    yCell = fileMPAS.variables["yCell"][:]
    zCell = fileMPAS.variables["zCell"][:]

    cellsOnCell = fileMPAS.variables["cellsOnCell"][:]
    nEdgesOnCell = fileMPAS.variables["nEdgesOnCell"][:]

    iceArea_CB_Initial   = fileMPAS.variables["iceAreaCell"][0,:]
    iceVolume_CB_Initial = fileMPAS.variables["iceVolumeCell"][0,:]
    iceThickness_CB_Initial = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_CB_Initial[iCell] > 0.0):
            iceThickness_CB_Initial[iCell] = iceVolume_CB_Initial[iCell] / iceArea_CB_Initial[iCell]

    iceArea_CB = fileMPAS.variables["iceAreaCell"][iTime,:]
    iceVolume_CB = fileMPAS.variables["iceVolumeCell"][iTime,:]
    iceThickness_CB = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_CB[iCell] > 0.0):
            iceThickness_CB[iCell] = iceVolume_CB[iCell] / iceArea_CB[iCell]

    fileMPAS.close()

    # slotted cylinder
    filename = "./output_%s_%s/output.2000.nc" %(experiment2,res)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    iceArea_SC_Initial   = fileMPAS.variables["iceAreaCell"][0,:]
    iceVolume_SC_Initial = fileMPAS.variables["iceVolumeCell"][0,:]
    iceThickness_SC_Initial = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_SC_Initial[iCell] > 0.0):
            iceThickness_SC_Initial[iCell] = iceVolume_SC_Initial[iCell] / iceArea_SC_Initial[iCell]

    iceArea_SC = fileMPAS.variables["iceAreaCell"][iTime,:]
    iceVolume_SC = fileMPAS.variables["iceVolumeCell"][iTime,:]
    iceThickness_SC = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_SC[iCell] > 0.0):
            iceThickness_SC[iCell] = iceVolume_SC[iCell] / iceArea_SC[iCell]

    print(np.amin(iceArea_SC_Initial), np.amax(iceArea_SC_Initial), np.amin(iceVolume_SC_Initial), np.amax(iceVolume_SC_Initial))
    print(np.amin(iceArea_SC), np.amax(iceArea_SC), np.amin(iceVolume_SC), np.amax(iceVolume_SC))

    fileMPAS.close()

    lons, yiceArea_CB_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB_Initial)
    lons, yiceArea_CB = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB)

    lons, yiceArea_SC_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC_Initial)
    lons, yiceArea_SC = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC)

    lons, yiceThickness_CB_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_CB_Initial)
    lons, yiceThickness_CB = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_CB)

    lons, yiceThickness_SC_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_SC_Initial)
    lons, yiceThickness_SC = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_SC)

    cm = 1/2.54  # centimeters in inches
    plt.rcParams["font.family"] = "Times New Roman"
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    linewidth = 1.0

    fig, axes = plt.subplots(1, 2, figsize=(15*cm,7*cm))

    # area
    axes[0].plot(lons, yiceArea_CB_Initial, linewidth=linewidth, color="black", ls="-")
    axes[0].plot(lons, yiceArea_CB,      linewidth=linewidth, color="red", ls="--")

    axes[0].legend(["Initial","MPM",""],frameon=False,fontsize=8)
    axes[0].set_xlabel("Longitude (radians)")
    axes[0].set_ylabel("Equatorial ice concentration")

    axes[0].set_xlim([0.5, math.pi-0.5])
    axes[0].set_title("(a) Cosine bell", loc='left')

    axes[1].plot(lons, yiceArea_SC_Initial, linewidth=linewidth, color="black", ls="-")
    axes[1].plot(lons, yiceArea_SC,      linewidth=linewidth, color="red", ls="--")

    axes[1].set_xlabel("Longitude (radians)")
    axes[1].set_xlim([0.5, math.pi-0.5])
    axes[1].set_title("(b) Slotted cylinder", loc='left')

    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("advection_equatorial.png",dpi=300)
    plt.savefig("advection_equatorial.eps")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    advection_equatorial()
