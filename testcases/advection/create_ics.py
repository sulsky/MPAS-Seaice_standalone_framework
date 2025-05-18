from netCDF4 import Dataset
import math, os, sys
import numpy as np
import numpy.ma as ma
from numba import njit

#--------------------------------------------------------------------

@njit
def slotted_cylinder(nCells,
                     xCell,
                     yCell,
                     zCell):

    iceAreaCell = np.zeros(nCells)
    iceVolumeCell = np.zeros(nCells)

    circleRadius = 0.5

    for iCell in range(0,nCells):

        r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

        if (r < circleRadius and yCell[iCell] > 0.0):

            iceAreaCell[iCell]   = 1.0
            iceVolumeCell[iCell] = 1.0

    for iCell in range(0,nCells):

        if (math.fabs(xCell[iCell]) < 1.0/12.0 and zCell[iCell] > -2.0/6.0):

            iceAreaCell[iCell]   = 0.0
            iceVolumeCell[iCell] = 0.0

    return iceAreaCell, iceVolumeCell

#--------------------------------------------------------------------

@njit
def cosine_bell_volume(nCells,
                       xCell,
                       yCell,
                       zCell):

    iceAreaCell = np.zeros(nCells)
    iceVolumeCell = np.zeros(nCells)

    circleRadius = 1.0/3.0

    for iCell in range(0,nCells):

        r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

        if (r < circleRadius and yCell[iCell] > 0.0):

            iceAreaCell[iCell]   = 1.0

            iceVolumeCell[iCell] = 0.5 * (1.0 + math.cos((math.pi * r) / circleRadius))

    return iceAreaCell, iceVolumeCell

#--------------------------------------------------------------------

@njit
def cosine_bell(nCells,
                xCell,
                yCell,
                zCell):

    iceAreaCell = np.zeros(nCells)
    iceVolumeCell = np.zeros(nCells)

    circleRadius = 1.0/3.0

    for iCell in range(0,nCells):

        r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

        if (r < circleRadius and yCell[iCell] > 0.0):

            iceAreaCell[iCell] = 0.5 * (1.0 + math.cos((math.pi * r) / circleRadius))

            iceVolumeCell[iCell] = 1.0

    return iceAreaCell, iceVolumeCell

#--------------------------------------------------------------------

@njit
def cylinder(nCells,
             xCell,
             yCell,
             zCell):

    iceAreaCell = np.zeros(nCells)
    iceVolumeCell = np.zeros(nCells)

    circleRadius = 0.5

    for iCell in range(0,nCells):

        r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

        if (r < circleRadius and yCell[iCell] > 0.0):

            iceAreaCell[iCell]   = 1.0
            iceVolumeCell[iCell] = 1.0

    return iceAreaCell, iceVolumeCell

#--------------------------------------------------------------------

def create_ic_file(res, icType, angle):

    # grid in
    gridFilename = "grid.%s.nc" %(res)
    gridFile = Dataset(gridFilename, "r")

    nCells = len(gridFile.dimensions["nCells"])
    nVertices = len(gridFile.dimensions["nVertices"])

    xCell = ma.getdata(gridFile.variables["xCell"][:])
    yCell = ma.getdata(gridFile.variables["yCell"][:])
    zCell = ma.getdata(gridFile.variables["zCell"][:])

    latVertex = ma.getdata(gridFile.variables["latVertex"][:])
    lonVertex = ma.getdata(gridFile.variables["lonVertex"][:])

    gridFile.close()

    # ice out
    icFilename = "ic_%s_%s.nc" %(icType, res)

    if (not os.path.isfile(icFilename)):

        icFile = Dataset(icFilename, "w", format="NETCDF3_CLASSIC")

        icFile.createDimension("nCells", nCells)
        icFile.createDimension("nVertices", nVertices)
        icFile.createDimension("nCategories", size=1)
        icFile.createDimension("ONE", size=1)

        uVelocity = icFile.createVariable("uVelocity", 'd', dimensions=("nVertices"))
        vVelocity = icFile.createVariable("vVelocity", 'd', dimensions=("nVertices"))

        days = 10.0
        seconds = days * 24.0 * 3600.0
        radius = 6371229.0
        uVelocityEquator = (2.0 * math.pi * radius) / (seconds)
        print("uVelocityEquator: ",uVelocityEquator)
        print("advection angle wrt equator: ",angle)
        for iVertex in range(0,nVertices):

            uVelocity[iVertex] =  uVelocityEquator * (math.cos(latVertex[iVertex]) * math.cos(angle)
                               + math.sin(latVertex[iVertex]) * math.cos(lonVertex[iVertex]) * math.sin(angle))
            vVelocity[iVertex] = -uVelocityEquator * (math.sin(lonVertex[iVertex]) * math.sin(angle))

        iceAreaCell   = icFile.createVariable("iceAreaCell",   'd', dimensions=("nCells"))
        iceVolumeCell = icFile.createVariable("iceVolumeCell", 'd', dimensions=("nCells"))

        iceAreaCategory   = icFile.createVariable("iceAreaCategory",   'd', dimensions=("nCells","nCategories","ONE"))
        iceVolumeCategory = icFile.createVariable("iceVolumeCategory", 'd', dimensions=("nCells","nCategories","ONE"))

        if (icType == "slotted_cylinder"):

            iceAreaCellArray, iceVolumeCellArray = slotted_cylinder(
                nCells,
                xCell,
                yCell,
                zCell)

        elif (icType == "cosine_bell_volume"):

            iceAreaCellArray, iceVolumeCellArray = cosine_bell_volume(
                nCells,
                xCell,
                yCell,
                zCell)

        elif (icType == "cosine_bell"):

            iceAreaCellArray, iceVolumeCellArray = cosine_bell(
                nCells,
                xCell,
                yCell,
                zCell)

        elif (icType == "cylinder"):

            iceAreaCellArray, iceVolumeCellArray = cylinder(
                nCells,
                xCell,
                yCell,
                zCell)

        iceAreaCell[:] = iceAreaCellArray[:]
        iceVolumeCell[:] = iceVolumeCellArray[:]

        iceAreaCategory[:,0]   = iceAreaCellArray[:]
        iceVolumeCategory[:,0] = iceVolumeCellArray[:]

        icFile.close()

#--------------------------------------------------------------------

def create_ics(angle = 0):

    reses = ["2562","10242","40962","163842"]

    icTypes = ["cosine_bell","slotted_cylinder"]

    for icType in icTypes:

        print("icType: ", icType)

        for res in reses:

            print("  Res: ", res)

            create_ic_file(res, icType, angle)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ics()
