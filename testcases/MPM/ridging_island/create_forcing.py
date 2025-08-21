from netCDF4 import Dataset
import netCDF4
from datetime import datetime, timedelta
from math import exp, sqrt, pow, sin, cos, radians
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import argparse

#-------------------------------------------------------------------------------

def create_forcing(gridFilename):

    # grid data
    fileGrid = Dataset(gridFilename,"r")

    nx = fileGrid.nx
    ny = fileGrid.ny

    nCells = len(fileGrid.dimensions["nCells"])

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:] - 1

    fileGrid.close()

    # atmos forcing
    nTimes = 120
    time0 = datetime(year=2000, month=1, day=1)
    time = time0
    dt = timedelta(hours=1)

    fileAtmos = Dataset("atmosphere_forcing_hourly.2000.nc","w",format="NETCDF3_CLASSIC")

    fileAtmos.createDimension("nCells", nCells)
    fileAtmos.createDimension("Time", None)
    fileAtmos.createDimension("StrLen", 64)

    uAirVelocity = np.zeros((nTimes,nCells))
    vAirVelocity = np.zeros((nTimes,nCells))
    xtimes = []

    for iTime in range(0, nTimes):

        timeStr = time.strftime("%Y-%m-%d_%H:%M:%S")
        elapsedTime = time - time0
        secs = elapsedTime.total_seconds()

        xtimes.append(timeStr)
        time = time + dt

        for iCell in range(0, nCells):

            x = xCell[iCell]
            y = yCell[iCell]

            uAirVelocity[iTime,iCell] = 10.0
            vAirVelocity[iTime,iCell] = 10.0

    varXtime = fileAtmos.createVariable("xtime","c",dimensions=["Time","StrLen"])
    uAirVelocityVar = fileAtmos.createVariable("uAirVelocity","d",dimensions=["Time","nCells"])
    vAirVelocityVar = fileAtmos.createVariable("vAirVelocity","d",dimensions=["Time","nCells"])

    for iTime in range(0,nTimes):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45
    uAirVelocityVar[:] = uAirVelocity[:]
    vAirVelocityVar[:] = vAirVelocity[:]

    fileAtmos.close()


    # ocean cyclone fields
    nTimes = 120
    time0 = datetime(year=2000, month=1, day=1)
    time = time0
    dt = timedelta(hours=1)

    fileOcean = Dataset("ocean_forcing_hourly.2000.nc","w",format="NETCDF3_CLASSIC")

    fileOcean.createDimension("nCells", nCells)
    fileOcean.createDimension("Time", None)
    fileOcean.createDimension("StrLen", 64)

    uOceanVelocity = np.zeros((nTimes,nCells))
    vOceanVelocity = np.zeros((nTimes,nCells))
    xtimes = []

    for iTime in range(0, nTimes):

        timeStr = time.strftime("%Y-%m-%d_%H:%M:%S")
        xtimes.append(timeStr)

        time = time + dt

        for iCell in range(0, nCells):

            x = xCell[iCell]
            y = yCell[iCell]

            uOceanVelocity[iTime,iCell] = 0.0
            vOceanVelocity[iTime,iCell] = 0.0

    varXtime = fileOcean.createVariable("xtime","c",dimensions=["Time","StrLen"])
    uOceanVelocityVar = fileOcean.createVariable("uOceanVelocity","d",dimensions=["Time","nCells"])
    vOceanVelocityVar = fileOcean.createVariable("vOceanVelocity","d",dimensions=["Time","nCells"])

    for iTime in range(0,nTimes):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45
    uOceanVelocityVar[:] = uOceanVelocity[:]
    vOceanVelocityVar[:] = vOceanVelocity[:]

    fileOcean.close()

    # plots
    wind = np.sqrt(np.add(np.power(uAirVelocity,2),
                          np.power(vAirVelocity,2)))

    currents = np.sqrt(np.add(np.power(uOceanVelocity,2),
                              np.power(vOceanVelocity,2)))

    nSkip = 8
    iTime = 48

    iCell = 0
    iCells = []
    for iCell in range(0,nCells):
        if (iCell % nSkip == 0):
            iCells.append(iCell)
        iCell += 1

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

    # winds
    fig, axis = plt.subplots()

    pc = PatchCollection(patches, cmap="jet")
    pc.set_array(np.array(wind[iTime,:]))
    axis.add_collection(pc)

    q = axis.quiver(xCell[iCells],
                    yCell[iCells],
                    uAirVelocity[iTime,iCells],
                    vAirVelocity[iTime,iCells],
                    wind[iTime,iCells])

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(pc, cax=cax)

    axis.set_aspect("equal")

    axis.set_xlim((xMin,xMax))
    axis.set_ylim((yMin,yMax))

    plt.tight_layout()
    plt.savefig("airVelocity.png",dpi=300)
    plt.cla()
    plt.close(fig)

    # currents
    fig, axis = plt.subplots()

    pc = PatchCollection(patches, cmap="jet")
    pc.set_array(np.array(currents[iTime,:]))
    axis.add_collection(pc)

    q = axis.quiver(xCell[iCells],
                    yCell[iCells],
                    uOceanVelocity[iTime,iCells],
                    vOceanVelocity[iTime,iCells],
                    currents[iTime,iCells])

    divider = make_axes_locatable(axis)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cb = fig.colorbar(q, cax=cax)

    axis.hlines(y=55000, xmin=0, xmax=100000, linewidth=0.5, color='r')
    axis.hlines(y=60000, xmin=0, xmax=100000, linewidth=0.5, color='r')
    axis.vlines(x=55000, ymin=0, ymax=100000, linewidth=0.5, color='r')
    axis.vlines(x=60000, ymin=0, ymax=100000, linewidth=0.5, color='r')

    axis.set_aspect("equal")

    axis.set_xlim((xMin,xMax))
    axis.set_ylim((yMin,yMax))

    plt.tight_layout()
    plt.savefig("oceanVelocity.png",dpi=300)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='gridFilename', required=True)

    args = parser.parse_args()

    create_forcing(args.gridFilename)
