from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import argparse

#-------------------------------------------------------------------------------

def create_ics(gridFilename):

    uAir = 10.0
    vAir =  0.0

    # grid data
    fileGrid = Dataset(gridFilename,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:] - 1
    cellsOnCell = fileGrid.variables["cellsOnCell"][:] - 1

    fileGrid.close()

    xmin = np.amin(xVertex)
    xmax = np.amax(xVertex)
    ymin = np.amin(yVertex)
    ymax = np.amax(yVertex)
    dx = xmax - xmin
    dy = ymax - ymin

    nCategories = 1

    iceAreaCategory = np.zeros((nCells,nCategories))
    iceVolumeCategory = np.zeros((nCells,nCategories))

    iceAreaCell = np.zeros((nCells))
    iceVolumeCell = np.zeros((nCells))

    iceThicknesses = np.zeros(nCategories)
    for iCategory in range(0,nCategories-1):
        iceThicknesses[iCategory] = 1.0
    #print(iceThicknesses)

    uAirVelocity = np.zeros(nCells)
    vAirVelocity = np.zeros(nCells)

    for iCell in range(0, nCells):

        x = xCell[iCell]
        y = yCell[iCell]
        uAirVelocity[iCell] = uAir
        vAirVelocity[iCell] = vAir

        if (x < 60000):
           uAirVelocity[iCell] = -uAir

        if (x < 100000 and x > 20000 and y < 100000 and y > 20000):
            for iCategory in range(0,nCategories):
                iceAreaCategory[iCell,iCategory] = 1.0
            iceAreaCell[iCell] = np.sum(iceAreaCategory[iCell,:])

            for iCategory in range(0,nCategories):
                iceVolumeCategory[iCell,iCategory] = iceThicknesses[iCategory] * iceAreaCategory[iCell,iCategory]
            iceVolumeCell[iCell] = np.sum(iceVolumeCategory[iCell,:])

    fileOut = Dataset("ic.nc", "w", format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("nCategories", nCategories)
    fileOut.createDimension("ONE", 1)

    var = fileOut.createVariable("iceAreaCategory", "d", dimensions=["nCells","nCategories","ONE"])
    var[:,:,0] = iceAreaCategory[:,:]

    var = fileOut.createVariable("iceVolumeCategory", "d", dimensions=["nCells","nCategories","ONE"])
    var[:,:,0] = iceVolumeCategory[:,:]

    var = fileOut.createVariable("iceAreaCell", "d", dimensions=["nCells"])
    var[:] = iceAreaCell[:]

    var = fileOut.createVariable("iceVolumeCell", "d", dimensions=["nCells"])
    var[:] = iceVolumeCell[:]

    fVertex = 1.46e-4
    var = fileOut.createVariable("fVertex", "d", dimensions=["nVertices"])
    var[:] = fVertex

    var = fileOut.createVariable("uAirVelocity", "d", dimensions=["nCells"])
    var[:] = uAirVelocity

    var = fileOut.createVariable("vAirVelocity", "d", dimensions=["nCells"])
    var[:] = vAirVelocity

    fileOut.close()

    # plots
    patches = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,4):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex],yVertex[iVertex]))
        patches.append(Polygon(vertices, closed=True, edgecolor="teal", fill=False, linewidth=0.1))

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    # iceAreaCell
    fig, axis = plt.subplots()

    pc = PatchCollection(patches, cmap="jet")
    pc.set_array(np.array(iceAreaCell))
    axis.add_collection(pc)

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pc, cax=cax)

    axis.set_aspect("equal")

    axis.set_xlim((xMin,xMax))
    axis.set_ylim((yMin,yMax))

    plt.tight_layout()
    plt.savefig("iceAreaCell.png",dpi=300)
    plt.cla()
    plt.close(fig)

    # iceVolumeCell
    fig, axis = plt.subplots()

    pc = PatchCollection(patches, cmap="jet")
    pc.set_array(np.array(iceVolumeCell))
    axis.add_collection(pc)

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pc, cax=cax)

    axis.set_aspect("equal")

    axis.set_xlim((xMin,xMax))
    axis.set_ylim((yMin,yMax))

    plt.tight_layout()
    plt.savefig("iceVolumeCell.png",dpi=300)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='gridFilename', required=True)

    args = parser.parse_args()

    create_ics(args.gridFilename)
