import sys

sys.path.append("../../../../utils/testcases")
from log_messages import log_message

import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def analyze_conservation():

    #----------------------------------------------
    # particles
    #----------------------------------------------
    filenames = sorted(glob.glob("output/part*"))

    totalParticleAreas   = []
    totalParticleVolumes = []
    for filename in filenames:

        filein = Dataset(filename,"r")

        nParticles = len(filein.dimensions["nParticles"])

        statusMP = filein.variables["statusMP"][0,:]

        areaMP = filein.variables["areaMP"][0,:]
        iceAreaCellMP = filein.variables["iceAreaCellMP"][0,:]
        iceVolumeCellMP = filein.variables["iceVolumeCellMP"][0,:]
        iceVolumeMP = filein.variables["iceVolumeCellMP"][0,:]

        filein.close()

        totalParticleArea   = 0.0
        totalParticleVolume = 0.0

        for iParticle in range(0,nParticles):
            if (statusMP[iParticle] == 1):

                totalParticleArea   += areaMP[iParticle] * iceAreaCellMP[iParticle]
                totalParticleVolume += areaMP[iParticle] * iceVolumeCellMP[iParticle]

        totalParticleAreas.append(totalParticleArea)
        totalParticleVolumes.append(totalParticleVolume)

    relParticleVolumeErrors = []
    for i in range(len(totalParticleAreas)):
        relParticleArea   = (totalParticleAreas  [i]-totalParticleAreas  [0])/totalParticleAreas  [0]
        relParticleVolume = (totalParticleVolumes[i]-totalParticleVolumes[0])/totalParticleVolumes[0]
        relParticleVolumeErrors.append(relParticleVolume)
        print("%5i: totalParticleArea: %.4f, rel: %10.7f, totalParticleVolume: %.4f, rel: %10.7f" \
              %(i, totalParticleAreas[i], relParticleArea, totalParticleVolumes[i], relParticleVolume))

    relParticleVolumeErrors = np.array(relParticleVolumeErrors)

    maxRelParticleVolumeError = np.amax(np.fabs(relParticleVolumeErrors))
    if (maxRelParticleVolumeError < 1e-12):
        message = "TEST PASSED: Max relative particle volume conservation error: %g" %(maxRelParticleVolumeError)
        log_message(message, "green", doPrint=True, logFile=None)
    else:
        message = "TEST WARNING: Max relative particle volume conservation error: %g" %(maxRelParticleVolumeError)
        log_message(message, "yellow", doPrint=True, logFile=None)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(totalParticleAreas)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Total ice area")
    axis.set_title("Total ice area on particles")

    plt.tight_layout()
    plt.savefig("particle_ice_area.png",dpi=300)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(totalParticleVolumes)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Total ice volume")
    axis.set_title("Total ice volume on particles")

    plt.tight_layout()
    plt.savefig("particle_ice_volume.png",dpi=300)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(relParticleVolumeErrors)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Relative error vs start")

    plt.tight_layout()
    plt.savefig("volume_conservation_particle.png",dpi=300)

    #----------------------------------------------
    # cells
    #----------------------------------------------
    filein = Dataset("./output/output.2000.nc","r")

    nTimes = len(filein.dimensions["Time"])
    nCells = len(filein.dimensions["nCells"])

    areaCell = filein.variables["areaCell"][:]
    iceAreaCell = filein.variables["iceAreaCell"][:]
    iceVolumeCell = filein.variables["iceVolumeCell"][:]

    filein.close()

    totalCellAreas = []
    totalCellVolumes = []
    for iTime in range(0,nTimes):
        totalCellArea = 0.0
        totalCellVolume = 0.0
        for iCell in range(0,nCells):
            totalCellArea   += areaCell[iCell] * iceAreaCell[iTime,iCell]
            totalCellVolume += areaCell[iCell] * iceVolumeCell[iTime,iCell]

        totalCellAreas.append(totalCellArea)
        totalCellVolumes.append(totalCellVolume)

    relCellVolumeErrors = []
    for i in range(len(totalCellAreas)):
        relCellArea   = (totalCellAreas  [i]-totalCellAreas  [0])/totalCellAreas  [0]
        relCellVolume = (totalCellVolumes[i]-totalCellVolumes[0])/totalCellVolumes[0]
        relCellVolumeErrors.append(relCellVolume)
        print("%5i: totalCellArea: %.4f, rel: %10.7f, totalCellVolume: %.4f, rel: %10.7f" \
              %(i, totalCellAreas[i], relCellArea, totalCellVolumes[i], relCellVolume))

    relCellVolumeErrors = np.array(relCellVolumeErrors)

    maxRelCellVolumeError = np.amax(np.fabs(relCellVolumeErrors))
    if (maxRelCellVolumeError < 1e-12):
        message = "TEST PASSED: Max relative cell volume conservation error: %g" %(maxRelCellVolumeError)
        log_message(message, "green", doPrint=True, logFile=None)
    else:
        message = "TEST FAIL: Max relative cell volume conservation error: %g" %(maxRelCellVolumeError)
        log_message(message, "red", doPrint=True, logFile=None)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(totalCellAreas)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Total ice area")
    axis.set_title("Total ice area on cells")

    plt.tight_layout()
    plt.savefig("cell_ice_area.png",dpi=300)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(totalCellVolumes)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Total ice volume")
    axis.set_title("Total ice volume on cells")

    plt.tight_layout()
    plt.savefig("cell_ice_volume.png",dpi=300)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(relCellVolumeErrors)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Relative error vs start")

    plt.tight_layout()
    plt.savefig("volume_conservation_cell.png",dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    analyze_conservation()
