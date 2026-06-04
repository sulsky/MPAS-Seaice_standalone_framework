from netCDF4 import Dataset, stringtochar
import numpy as np
import os
import argparse
from math import fabs
from tqdm import tqdm

#-------------------------------------------------------------------------------

def latlon_to_xyz(lat, lon, radius):

    x = radius * np.cos(lat) * np.cos(lon)
    y = radius * np.cos(lat) * np.sin(lon)
    z = radius * np.sin(lat)

    return x, y, z

#-------------------------------------------------------------------------------

def create_mpas_greenland_calving_rates(meshFilename,
                                        calvingFilename,
                                        regionType):

    # get input calving and region files
    MPAS_SEAICE_STANDALONE_DATA = os.environ.get('MPAS_SEAICE_STANDALONE_DATA')
    if (MPAS_SEAICE_STANDALONE_DATA is None):
        raise Exception("MPAS_SEAICE_STANDALONE_DATA must be set")

    filenameGate = "%s/icebergs/Greenland/gate.nc" %(MPAS_SEAICE_STANDALONE_DATA)
    if (not os.path.isfile(filenameGate)):
        raise Exception("Could not find calving file: %s" %(filenameGate))

    earthRadius = 6371000.0

    # load mpas mesh file
    fileMesh = Dataset(meshFilename,"r")

    nCells = len(fileMesh.dimensions["nCells"])

    latCell = fileMesh.variables["latCell"][:]
    lonCell = fileMesh.variables["lonCell"][:]

    xCell, yCell, zCell = latlon_to_xyz(latCell, lonCell, earthRadius)

    fileMesh.close()

    # gates
    fileGate = Dataset(filenameGate,"r")

    nGates = len(fileGate.dimensions["gate"])

    discharge = fileGate.variables["discharge"][:] # already Gt/y
    discharge = discharge.astype(np.double)

    lonGate = fileGate.variables["mean_lon"][:]
    latGate = fileGate.variables["mean_lat"][:]

    lonGate = np.radians(lonGate)
    latGate = np.radians(latGate)

    gateNames = fileGate.variables["name_Mouginot"][:]

    regionsMouginotByGate = fileGate.variables["region"][:]
    regionsMouginot = set()
    for iRegion in range(0,regionsMouginotByGate.shape[0]):
        regionsMouginot.add(regionsMouginotByGate[iRegion])
    regionsMouginot = list(regionsMouginot)
    nRegionsMouginot = len(regionsMouginot)
    regionIndexMap = {}
    for iRegion in range(0,nRegionsMouginot):
        regionIndexMap[regionsMouginot[iRegion]] = iRegion
    gateIndexMouginot = {}
    for iGate in range(0,nGates):
        gateIndexMouginot[regionsMouginotByGate[iGate]] = iGate

    fileGate.close()

    # average dischange over time
    dischargeAvg = np.mean(discharge,axis=1)
    dischargeTotal = np.sum(dischargeAvg)
    dischargeMouginot = np.zeros(nRegionsMouginot)
    for iGate in range(0,nGates):
        iRegion = regionIndexMap[regionsMouginotByGate[iGate]]
        dischargeMouginot[iRegion] += dischargeAvg[iGate]

    dischargeMouginotSum = np.sum(dischargeMouginot)
    if (fabs(dischargeMouginotSum-dischargeTotal) > 1e-11):
        raise Exception("Discharge incorrect 1 %f %f %f" %(dischargeMouginotSum,dischargeTotal,fabs(dischargeMouginotSum-dischargeTotal)))

    if (regionType == "gate"):
        nCalvingRegions = nGates
        calvingRegionNames = [str(x) for x in gateNames]
        calvingRateRegions = dischargeAvg
    elif (regionType == "Mouginot_Region"):
        nCalvingRegions = nRegionsMouginot
        calvingRegionNames = [str(x) for x in regionsMouginot]
        calvingRateRegions = dischargeMouginot
    elif (regionType is None):
        nCalvingRegions = 1
        calvingRegionNames = ["NONE"]
        calvingRateRegions = np.array([np.sum(dischargeAvg)])

    # find gate cell
    iCellMin = np.zeros(nGates,dtype="i")
    for iGate in tqdm(range(0,nGates)):
        xGate, yGate, zGate = latlon_to_xyz(latGate[iGate], lonGate[iGate], earthRadius)

        dist = np.add(np.add(np.power(xCell-xGate,2),
                             np.power(yCell-yGate,2)),
                             np.power(zCell-zGate,2))

        iCellMin[iGate] = np.argmin(dist)

    # max gates per cell
    regionIndices = [set() for _ in range(nCells)]
    for iGate in tqdm(range(0,nGates)):
        iCell = iCellMin[iGate]
        if (regionType == "gate"):
            regionIndices[iCell].add(iGate)
        elif (regionType == "Mouginot_Region"):
            print(iGate,regionsMouginotByGate[iGate],regionIndexMap[regionsMouginotByGate[iGate]])
            regionIndices[iCell].add(regionIndexMap[regionsMouginotByGate[iGate]])
        elif (regionType is None):
            regionIndices[iCell].add(0)

    nCalvingRegionsPerCell = np.zeros(nCells,dtype="i")
    for iCell in range(0,nCells):
        regionIndices[iCell] = list(regionIndices[iCell])
        nCalvingRegionsPerCell[iCell] = len(regionIndices[iCell])

    maxCalvingRegionsPerCell = np.amax(nCalvingRegionsPerCell)

    print("maxCalvingRegionsPerCell: ", maxCalvingRegionsPerCell)

    # output - Gt/y
    calvingRate = np.zeros((nCells,maxCalvingRegionsPerCell))
    for iGate in range(0,nGates):
        iCell = iCellMin[iGate]
        if (regionType == "gate"):
            regionIndex = iGate
        elif (regionType == "Mouginot_Region"):
            regionIndex = regionIndexMap[regionsMouginotByGate[iGate]]
        elif (regionType is None):
            regionIndex = 0
        regionIndexInCell = regionIndices[iCell].index(regionIndex)
        calvingRate[iCell,regionIndexInCell] += dischargeAvg[iGate]

    calvingRateSum = np.sum(calvingRate)
    if (fabs(calvingRateSum - dischargeTotal) > 1e-11):
        raise Exception("Discharge incorrect 2 %f %f %f" %(calvingRateSum,dischargeTotal,fabs(calvingRateSum-dischargeTotal)))

    # gate index
    calvingRegionIndex = np.ones((nCells,maxCalvingRegionsPerCell),dtype="i")
    calvingRegionIndex[:] *= -1
    for iCell in range(0,nCells):
        regionIndex = list(regionIndices[iCell])
        for iRegionPerCell in range(0,nCalvingRegionsPerCell[iCell]):
            calvingRegionIndex[iCell,iRegionPerCell] = regionIndex[iRegionPerCell]

    # output file
    fileOut = Dataset(calvingFilename,"w",format="NETCDF3_CLASSIC")

    if (regionType == "gate"):
        fileOut.regionType = "gate"
    elif (regionType == "Mouginot_Region"):
        fileOut.regionType = "Mouginot_Region"
    elif (regionType is None):
        fileOut.regionType = "none"

    fileOut.totalCalvingRate = dischargeTotal
    fileOut.src = "doi:10.5194/essd-12-1367-2020 doi:10.22008/promice/data/ice_discharge"

    fileOut.createDimension("nCells",nCells)
    fileOut.createDimension("maxCalvingRegionsPerCell",maxCalvingRegionsPerCell)
    fileOut.createDimension("nCalvingRegions",nCalvingRegions)
    StrLen = 64
    fileOut.createDimension("StrLen",StrLen)

    var = fileOut.createVariable("nCalvingRegionsPerCell","i",dimensions=["nCells"])
    var[:] = nCalvingRegionsPerCell[:]

    var = fileOut.createVariable("calvingRate","d",dimensions=["nCells","maxCalvingRegionsPerCell"])
    var.units = "Gt/y"
    var[:] = calvingRate[:]

    var = fileOut.createVariable("calvingRegionIndex","i",dimensions=["nCells","maxCalvingRegionsPerCell"])
    var[:] = calvingRegionIndex[:]

    var = fileOut.createVariable("calvingRegionNames","S1",dimensions=["nCalvingRegions","StrLen"])
    char_array = np.array(
        [list(s.ljust(StrLen)) for s in calvingRegionNames],
        dtype="S1"
    )
    var[:] = char_array

    var = fileOut.createVariable("calvingRateRegions", "d", dimensions=["nCalvingRegions"])
    var.units = "Gt/y"
    var[:] = calvingRateRegions[:]

    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m', dest='meshFilename', required=True, help='Create MPAS-Seaice input iceberg calving file for Greenland')
    parser.add_argument('-o', dest='calvingFilename', default="calving_mpas_greenland.nc", help='MPAS calving input file name')
    parser.add_argument('-r', dest='regionType', default=None, choices=["gate","Mouginot_Region"], help='List calving by source calving region')

    args = parser.parse_args()

    create_mpas_greenland_calving_rates(args.meshFilename,
                                        args.calvingFilename,
                                        args.regionType)
