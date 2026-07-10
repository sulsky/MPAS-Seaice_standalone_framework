import sys

sys.path.append("../../../../utils/testcases")
from log_messages import log_message

import glob
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def analyze_conservation():

    filenames = sorted(glob.glob("output/part*"))

    totalAreas   = []
    totalVolumes = []
    for filename in filenames:

        filein = Dataset(filename,"r")

        nParticles = len(filein.dimensions["nParticles"])

        statusMP = filein.variables["statusMP"][0,:]

        areaMP = filein.variables["areaMP"][0,:]
        iceAreaCellMP = filein.variables["iceAreaCellMP"][0,:]
        iceVolumeCellMP = filein.variables["iceVolumeCellMP"][0,:]
        iceVolumeMP = filein.variables["iceVolumeCellMP"][0,:]

        filein.close()

        totalArea   = 0.0
        totalVolume = 0.0

        for iParticle in range(0,nParticles):
            if (statusMP[iParticle] == 1):

                totalArea   += areaMP[iParticle] * iceAreaCellMP[iParticle]
                totalVolume += areaMP[iParticle] * iceVolumeCellMP[iParticle]

        print(filename,":, totalArea: ", totalArea, ", totalVolume: ", totalVolume)

        totalAreas.append(totalArea)
        totalVolumes.append(totalVolume)

    relVolumeErrors = []

    for i in range(len(totalAreas)):
        relVolumeError = (totalVolumes[i]-totalVolumes[0])/totalVolumes[0]
        relVolumeErrors.append(relVolumeError)

    relVolumeErrors = np.array(relVolumeErrors)

    fig, axis = plt.subplots()

    axis.axhline(y=0, color='grey', linewidth=1)
    axis.plot(relVolumeErrors)

    axis.set_xlabel("Output index")
    axis.set_ylabel("Relative error vs start")

    plt.savefig("volume_conservation.png",dpi=300)

    maxRelVolumeError = np.amax(relVolumeErrors)

    if (maxRelVolumeError < 1e-12):
        message = "TEST PASSED: Max relative volume conservation error: %g" %(maxRelVolumeError)
        log_message(message, "green", logFile=None)
    else:
        message = "TEST FAILED: Max relative volume conservation error: %g" %(maxRelVolumeError)
        log_message(message, "red", logFile=None)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    analyze_conservation()
