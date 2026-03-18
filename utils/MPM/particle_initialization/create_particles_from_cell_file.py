from initial_particle_positions import place_particles_in_cell, create_particles_file, plot_particles

from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import argparse

#-------------------------------------------------------------------------------

def create_particles_from_cell_file(filenameMesh,
                                    filenameSeaice,
                                    particleInitType,
                                    particleInitNumber,
                                    posnInitType,
                                    filenameOut):

    # mesh info
    fileMesh = Dataset(filenameMesh,"r")

    on_a_sphere = fileMesh.on_a_sphere
    sphere = fileMesh.on_a_sphere
    if (sphere == "NO"):
        on_a_sphere = False
    elif (sphere == "YES"):
        on_a_sphere = True

    nCells = len(fileMesh.dimensions["nCells"])

    nEdgesOnCell   = ma.getdata(fileMesh.variables["nEdgesOnCell"][:])
    verticesOnCell = ma.getdata(fileMesh.variables["verticesOnCell"][:])
    indexToCellID  = ma.getdata(fileMesh.variables["indexToCellID"][:])

    areaCell  = ma.getdata(fileMesh.variables["areaCell"][:])
    latVertex = ma.getdata(fileMesh.variables["latVertex"][:])
    lonVertex = ma.getdata(fileMesh.variables["lonVertex"][:])
    xVertex   = ma.getdata(fileMesh.variables["xVertex"][:])
    yVertex   = ma.getdata(fileMesh.variables["yVertex"][:])
    zVertex   = ma.getdata(fileMesh.variables["zVertex"][:])
    xCell     = ma.getdata(fileMesh.variables["xCell"][:])
    yCell     = ma.getdata(fileMesh.variables["yCell"][:])
    zCell     = ma.getdata(fileMesh.variables["zCell"][:])

    fileMesh.close()

    verticesOnCell[:] -= 1

    # average cell size
    averageCellSize = np.mean(areaCell)

    # particle positions
    nParticlesCell = np.zeros(nCells,dtype="i")

    # Eulerian ice area cell
    fileIn = Dataset(filenameSeaice,"r")

    nCategories = len(fileIn.dimensions["nCategories"])

    iceAreaCell       = fileIn.variables["iceAreaCell"][:]
    iceVolumeCell     = fileIn.variables["iceVolumeCell"][:]
    iceAreaCategory   = fileIn.variables["iceAreaCategory"][:]
    iceVolumeCategory = fileIn.variables["iceVolumeCategory"][:]

    fileIn.close()

    posnMP = []
    latCellMP = []
    lonCellMP = []
    areaMP = []
    cellIDCreationMP = []
    creationIndexMP = []
    iceAreaCellMP = []
    iceVolumeCellMP = []
    iceAreaCategoryMP = []
    iceVolumeCategoryMP = []

    for iCell in range(0, nCells):

        if (iceAreaCell[iCell] > 0.0):

            xInCell, \
                yInCell, \
                zInCell, \
                areaInCell = place_particles_in_cell(particleInitType,
                                                  particleInitNumber,
                                                  areaCell[iCell],
                                                  averageCellSize,
                                                  posnInitType,
                                                  on_a_sphere,
                                                  0.0,
                                                  nEdgesOnCell[iCell],
                                                  verticesOnCell[iCell,:],
                                                  latVertex,
                                                  lonVertex,
                                                  xVertex,
                                                  yVertex,
                                                  zVertex,
                                                  xCell[iCell],
                                                  yCell[iCell],
                                                  zCell[iCell])

            k = 0
            for x, y, z, area in zip(xInCell, yInCell, zInCell, areaInCell):
                posnMP.append([x, y, z])
                latCellMP.append(0.0)
                lonCellMP.append(0.0)
                areaMP.append(area)
                cellIDCreationMP.append(indexToCellID[iCell])
                creationIndexMP.append(k)
                iceAreaCellMP.append(iceAreaCell[iCell])
                iceVolumeCellMP.append(iceVolumeCell[iCell])
                iceAreaCategoryMP.append(iceAreaCategory[iCell])
                iceVolumeCategoryMP.append(iceVolumeCategory[iCell])
                k = k + 1

    nParticles = len(posnMP)
    posnMP = np.array(posnMP)
    latCellMP = np.array(latCellMP)
    lonCellMP = np.array(lonCellMP)
    areaMP = np.array(areaMP)
    cellIDCreationMP = np.array(cellIDCreationMP)
    creationIndexMP = np.array(creationIndexMP)
    nParticlesCell = np.array(nParticlesCell)
    iceAreaCellMP = np.array(iceAreaCellMP)
    iceVolumeCellMP = np.array(iceVolumeCellMP)
    iceAreaCategoryMP = np.array(iceAreaCategoryMP)
    iceVolumeCategoryMP = np.array(iceVolumeCategoryMP)

    # output
    create_particles_file(filenameOut,
                          nParticles,
                          nCells,
                          nCategories,
                          posnMP,
                          latCellMP,
                          lonCellMP,
                          areaMP,
                          cellIDCreationMP,
                          creationIndexMP,
                          nParticlesCell,
                          iceAreaCellMP,
                          iceVolumeCellMP,
                          iceAreaCategoryMP,
                          iceVolumeCategoryMP)

    # plot
    plot_particles(on_a_sphere,
                   posnMP,
                   xVertex,
                   yVertex,
                   zVertex,
                   0.0,
                   "particles.png")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest='filenameMesh', required=True)
    parser.add_argument('-i', dest='filenameSeaice', required=True)
    parser.add_argument('-t', dest="particleInitType", required=True, choices=["number","area"])
    parser.add_argument('-n', dest="particleInitNumber", required=True, type=int)
    parser.add_argument('-p', dest="particlePositionInitType", required=True, choices=["even","onePerEdge","poisson","random"])
    parser.add_argument('-o', dest='filenameOut', required=True)

    args = parser.parse_args()

    create_particles_from_cell_file(args.filenameMesh,
                                    args.filenameSeaice,
                                    args.particleInitType,
                                    args.particleInitNumber,
                                    args.particlePositionInitType,
                                    args.filenameOut)
