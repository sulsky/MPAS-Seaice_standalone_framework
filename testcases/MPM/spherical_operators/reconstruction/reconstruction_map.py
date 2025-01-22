from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math

degreesToRadians = math.pi / 180.0

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#---------------------------------------------------------------

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > 20.0*degreesToRadians):

            polygonVertices = []

            useVertex = False
            for iCellOnVertex in range(0,vertexDegree):

                iCell = cellsOnVertex[iVertex,iCellOnVertex] - 1

                polygonVertices.append((xCell[iCell],yCell[iCell]))

                if (xCell[iCell] >= minX and xCell[iCell] <= maxX and \
                    yCell[iCell] >= minY and yCell[iCell] <= maxY):
                    useVertex = True

            if (useVertex):
                polygon = Polygon(polygonVertices)
                patches.append(polygon)

                colours.append(mpasArray[iVertex])

                minval = min(minval,mpasArray[iVertex])
                maxval = max(maxval,mpasArray[iVertex])

    patchCollection = PatchCollection(patches, cmap=cmap, rasterized=True)
    patchCollection.set_array(np.array(colours))
    patchCollection.set_linewidth(0)

    patchCollection.set_clim(vmin=vmin,vmax=vmax)

    return patchCollection, minval, maxval

#---------------------------------------------------------------

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
    axes.add_collection(patchCollection)
    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        axes.set_title(title, fontsize=8)

    if (subfigureLabel != None):
        axes.text(0.02, 0.89, subfigureLabel, verticalalignment='bottom', horizontalalignment='left',transform=axes.transAxes, fontsize=8)

    if (colorbar):
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = fig.colorbar(patchCollection,cax=cax)
        if (unityBar):
            cb.ax.set_yticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
        if (sciNote):
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()

#---------------------------------------------------------------

def reconstruction_map():

    # grid
    fileGrid = Dataset("grid.40962.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    # ic
    fileIC = Dataset("ic_latlon_40962.nc","r")

    uVelocityAnalytical = fileIC.variables["uVelocity"][:]
    vVelocityAnalytical = fileIC.variables["vVelocity"][:]

    print("Velocity: ",
          np.amin(uVelocityAnalytical), np.amax(uVelocityAnalytical),
          np.amin(vVelocityAnalytical), np.amax(vVelocityAnalytical))

    fileIC.close()

    # mpm
    fileMPM = Dataset("./output_latlon_40962/output.2000.nc","r")

    uVelocityMPM = fileMPM.variables["uVelocity"][0,:]
    vVelocityMPM = fileMPM.variables["vVelocity"][0,:]

    uVelocityDiff = (uVelocityMPM - uVelocityAnalytical)
    vVelocityDiff = (vVelocityMPM - vVelocityAnalytical)

    print("MPM Interpolation Error:  ",
          np.amin(uVelocityDiff), np.amax(uVelocityDiff),
          np.amin(vVelocityDiff), np.amax(vVelocityDiff))

    fileMPM.close()


    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    minVelDiff = -1.0e-7
    maxVelDiff =  1.0e-7

    fig, axes = plt.subplots(2, 2)

    fig.set_size_inches(7, 6.75)

    plot_subfigure(axes[0,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocityAnalytical, minVelocity, maxVelocity, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocityAnalytical, minVelocity, maxVelocity, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(b)', False)

    plot_subfigure(axes[1,0], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocityDiff, minVelDiff, maxVelDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM diff ($u^\prime$ direction)', '(c)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocityDiff, minVelDiff, maxVelDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM diff ($v^\prime$ direction)', '(d)', False)

    #plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("reconstruction_map.png",dpi=400)

    plt.clf()
    plt.cla()
    plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    reconstruction_map()
