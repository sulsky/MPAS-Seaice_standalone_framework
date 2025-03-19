from netCDF4 import Dataset
import numpy as np

#---------------------------------------------------------

def average_particle_ice_area_to_cell(filenameParticles, nCells):

    filein = Dataset(filenameParticles,"r")

    nParticles = len(filein.dimensions["nParticles"])

    statusMP      = filein.variables["statusMP"][0,:]
    iCellMP       = filein.variables["iCellMP"][0,:]-1
    iceAreaCellMP = filein.variables["iceAreaCellMP"][0,:]

    filein.close()

    iceAreaCell = np.zeros(nCells)
    nParticlesCell = np.zeros(nCells)

    for iParticle in range(0,nParticles):
        if (statusMP[iParticle] == 1):
            iCell = iCellMP[iParticle]
            iceAreaCell[iCell] += iceAreaCellMP[iParticle]
            nParticlesCell[iCell] += 1.0

    for iCell in range(0,nCells):
        if (nParticlesCell[iCell] > 0.0):
            iceAreaCell[iCell] /= nParticlesCell[iCell]

    return iceAreaCell

#---------------------------------------------------------
