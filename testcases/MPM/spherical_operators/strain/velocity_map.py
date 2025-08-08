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
            for iCellOnVertex in range(0,vertexDegree[iVertex]):

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

def velocity_map():

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
    fileIC = Dataset("ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    print("Velocity: ",
          np.amin(uVelocity), np.amax(uVelocity),
          np.amin(vVelocity), np.amax(vVelocity))

    # mpm
    fileMPM = Dataset("./output_mpmvar_40962/output.2000.nc","r")

    uVelocityMPM = fileMPM.variables["uVelocityInitial"][0,:]
    vVelocityMPM = fileMPM.variables["vVelocityInitial"][0,:]

    uVelocityDiff = uVelocity - uVelocityMPM
    vVelocityDiff = vVelocity - vVelocityMPM

    print("MPM Diff: ",
          np.amin(uVelocityDiff), np.amax(uVelocityDiff),
          np.amin(vVelocityDiff), np.amax(vVelocityDiff))

    fileMPM.close()

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(3, 2)
    fig.set_size_inches(6, 4)

    minVelo = -3.3
    maxVelo =  3.3

    minDiff = -5.e-4
    maxDiff =  5.e-4

    # Velocities
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(b)', True)

    # MPM
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocityDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$ MPM.', '(c)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocityDiff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$ MPM.', '(d)', True)
    plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocityDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$u^\prime$ MPM.', '(e)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocityDiff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$v^\prime$ MPM.', '(f)', True)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("velocity_map.png",dpi=400)
    #plt.savefig("velocity_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()

#---------------------------------------------------------------

if __name__ == "__main__":

    velocity_map()
