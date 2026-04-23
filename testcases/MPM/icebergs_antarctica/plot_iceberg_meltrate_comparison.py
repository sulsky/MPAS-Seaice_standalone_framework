import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.colors as mcolors
import glob
import calendar
import argparse
import os
import sys

#-------------------------------------------------------------------------------

def meltrate_comparison(month):

    if (month == 0):
        dateStr = "All year"
    else:
        dateStr = calendar.month_abbr[month]

    # load mesh data
    fileMesh = Dataset("grid.nc", "r")

    nCells = len(fileMesh.dimensions["nCells"])
    nEdges = len(fileMesh.dimensions["nEdges"])

    latCell = fileMesh.variables["latCell"][:]
    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]-1
    cellsOnEdge = fileMesh.variables["cellsOnEdge"][:]-1
    verticesOnEdge = fileMesh.variables["verticesOnEdge"][:]-1
    latCell = fileMesh.variables["latCell"][:]
    lonCell = fileMesh.variables["lonCell"][:]
    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]
    latEdge = fileMesh.variables["latEdge"][:]
    lonEdge = fileMesh.variables["lonEdge"][:]
    areaCell = fileMesh.variables["areaCell"][:]

    fileMesh.close()

    latVertex = np.degrees(latVertex)
    lonVertex = np.degrees(lonVertex)
    latCell = np.degrees(latCell)
    lonCell = np.degrees(lonCell)
    for iCell in range(0,nCells):
        if (lonCell[iCell] > 180.0):
            lonCell[iCell] -= 360.0
    latEdge = np.degrees(latEdge)
    lonEdge = np.degrees(lonEdge)
    for iEdge in range(0,nEdges):
        if (lonEdge[iEdge] > 180.0):
            lonEdge[iEdge] -= 360.0

    boundaryEdge = np.zeros(nEdges,dtype="i")
    for iEdge in range(0,nEdges):
        if (cellsOnEdge[iEdge,0] == -1 or
            cellsOnEdge[iEdge,1] == -1):
            boundaryEdge[iEdge] = 1

    lineSegments = []
    for iEdge in range(0,nEdges):
        if (boundaryEdge[iEdge] == 1):
            iVertex1 = verticesOnEdge[iEdge,0]
            iVertex2 = verticesOnEdge[iEdge,1]
            lonVertexUse1 = lonVertex[iVertex1]
            if (lonEdge[iEdge] < -90.0 and
                lonVertex[iVertex1] > 90.0):
                lonVertexUse1 = lonVertexUse1 - 360.0
            elif (lonEdge[iEdge] > 90.0 and
                  lonVertex[iVertex1] < -90.0):
                lonVertexUse1 = lonVertexUse1 + 360.0
            lonVertexUse2 = lonVertex[iVertex2]
            if (lonEdge[iEdge] < -90.0 and
                lonVertex[iVertex2] > 90.0):
                lonVertexUse2 = lonVertexUse2 - 360.0
            elif (lonEdge[iEdge] > 90.0 and
                  lonVertex[iVertex2] < -90.0):
                lonVertexUse2 = lonVertexUse2 + 360.0
            lineSegments.append([[lonVertexUse1,latVertex[iVertex1]],
                                 [lonVertexUse2,latVertex[iVertex2]]])

    lc = LineCollection(lineSegments, color="black", linestyle='solid', linewidth=0.2)


    icebergMeltRateCellAvg = np.zeros(nCells)
    nAvg = 0

    if (month == 0):
        filenameMPASTemplates = "analysis_members/timeSeriesStatsMonthly.*.nc"
    else:
        filenameMPASTemplates = "analysis_members/timeSeriesStatsMonthly.*-%2.2i.nc" %(month)
    filenamesMPAS = sorted(glob.glob(filenameMPASTemplates))
    for filenameMPAS in filenamesMPAS:

        fileMPAS = Dataset(filenameMPAS,"r")

        icebergMeltRateCell = fileMPAS.variables["timeMonthly_avg_icebergMeltRateCell"][0,:]

        icebergDensity = fileMPAS.config_iceberg_density

        fileMPAS.close()

        icebergMeltRateCellAvg[:] += icebergMeltRateCell[:]
        nAvg += 1

    icebergMeltRateCellAvg[:] /= float(nAvg)

    icebergMeltRateCellAvg[:] /= areaCell[:]
    icebergMeltRateCellAvg[:] *= icebergDensity

    vminMPAS =  sys.float_info.max
    vmaxMPAS = -sys.float_info.max

    patchesMPAS = []
    colorsMPAS = []
    for iCell in range(0,nCells):
        if (icebergMeltRateCellAvg[iCell] > 0.0):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                lonVertexUse = lonVertex[iVertex]
                if (lonCell[iCell] < -90.0 and
                    lonVertex[iVertex] > 90.0):
                    lonVertexUse = lonVertexUse - 360.0
                elif (lonCell[iCell] > 90.0 and
                    lonVertex[iVertex] < -90.0):
                    lonVertexUse = lonVertexUse + 360.0
                vertices.append([lonVertexUse,latVertex[iVertex]])
            patchesMPAS.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))
            colorsMPAS.append(icebergMeltRateCellAvg[iCell])
            vminMPAS = min(vminMPAS,icebergMeltRateCellAvg[iCell])
            vmaxMPAS = min(vmaxMPAS,icebergMeltRateCellAvg[iCell])

    # Merino data
    MPAS_SEAICE_STANDALONE_DATA = os.environ.get('MPAS_SEAICE_STANDALONE_DATA')
    if (MPAS_SEAICE_STANDALONE_DATA is None):
        raise Exception("MPAS_SEAICE_STANDALONE_DATA must be set")

    filein = Dataset("%s/icebergs/mmc2.nc" %(MPAS_SEAICE_STANDALONE_DATA),"r")

    nx = len(filein.dimensions["x"])
    ny = len(filein.dimensions["y"])

    longitude = filein.variables["longitude"][:]
    latitude = filein.variables["latitude"][:]
    Icb_flux = filein.variables["Icb_flux"][:]

    filein.close()

    longitudeMerino = np.concatenate((longitude[0,430:],longitude[0,0:430]))
    latitudeMerino = latitude[:,0]
    icebergMeltfluxMerino = np.concatenate((Icb_flux[:,:,430:],Icb_flux[:,:,0:430]),axis=2)

    if (month == 0):
        icebergMeltfluxMerino = np.mean(icebergMeltfluxMerino,axis=0)
    else:
        icebergMeltfluxMerino = icebergMeltfluxMerino[month-1,:,:]

    vminMerino = np.min(icebergMeltfluxMerino[np.nonzero(icebergMeltfluxMerino)])
    vmaxMerino = np.max(icebergMeltfluxMerino[np.nonzero(icebergMeltfluxMerino)])

    vmin = min(vminMPAS,vminMerino)
    vmax = max(vmaxMPAS,vmaxMerino)

    # plot
    plt.rcParams["font.family"] = "Times New Roman"

    fig, axes = plt.subplots(2,1)

    # MPAS plot
    pcMPAS = PatchCollection(patchesMPAS, match_original=True, cmap="jet", norm=mcolors.LogNorm(vmin=vmin, vmax=vmax))
    pcMPAS.set_array(colorsMPAS)
    axes[0].add_collection(pcMPAS)
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes('right', size='2%', pad=0.02)
    cb = fig.colorbar(pcMPAS, cax=cax)
    cb.set_label("kg/m2/s")

    lcMPAS = LineCollection(lineSegments, color="black", linestyle='solid', linewidth=0.2)
    axes[0].add_collection(lcMPAS)

    axes[0].autoscale_view()

    axes[0].set_xlim(-180.0,180.0)
    axes[0].set_ylim(-80.0,-40.0)

    axes[0].set_title("MPAS-Seaice - %s" %(dateStr))
    
    # Merino plot
    sc = axes[1].pcolormesh(longitudeMerino,
                            latitudeMerino,
                            icebergMeltfluxMerino,
                            shading='nearest', cmap='jet', norm=mcolors.LogNorm(vmin=vmin, vmax=vmax))
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes('right', size='2%', pad=0.02)
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("kg/m2/s")

    lcMerino = LineCollection(lineSegments, color="black", linestyle='solid', linewidth=0.2)
    axes[1].add_collection(lcMerino)

    axes[1].autoscale_view()

    axes[1].set_xlim(-180.0,180.0)
    axes[1].set_ylim(-80.0,-40.0)

    axes[1].set_title("Merino - %s" %(dateStr))

    plt.tight_layout()
    plt.savefig("meltrate.png",dpi=600)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m', dest='month', default=0, type=int, help='')

    args = parser.parse_args()

    meltrate_comparison(args.month)
