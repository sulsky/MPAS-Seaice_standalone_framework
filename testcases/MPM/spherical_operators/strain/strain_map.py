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

def strain_map():

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

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    print("Strain: ",
          np.amin(strain11VertexAnalytical), np.amax(strain11VertexAnalytical),
          np.amin(strain22VertexAnalytical), np.amax(strain22VertexAnalytical),
          np.amin(strain12VertexAnalytical), np.amax(strain12VertexAnalytical))

    fileIC.close()

    # mpm
    fileMPM = Dataset("./output_mpm_40962/output.2000.nc","r")

    strain11 = fileMPM.variables["strain11AvgVertex"][0,:]
    strain22 = fileMPM.variables["strain22AvgVertex"][0,:]
    strain12 = fileMPM.variables["strain12AvgVertex"][0,:]

    strain11Diff = strain11 - strain11VertexAnalytical
    strain22Diff = strain22 - strain22VertexAnalytical
    strain12Diff = strain12 - strain12VertexAnalytical

    print("MPM: ",
          np.amin(strain11Diff), np.amax(strain11Diff),
          np.amin(strain22Diff), np.amax(strain22Diff),
          np.amin(strain12Diff), np.amax(strain12Diff))

    fileMPM.close()

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(6, 6)
    fig.set_size_inches(9, 8)

    minStrain = -3.3
    maxStrain =  3.3

    minDiff = -0.043
    maxDiff =  0.043

    # Velocities
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(b)', True)
    axes[0,2].axis('off')
    axes[0,3].axis('off')
    axes[0,4].axis('off')
    axes[0,5].axis('off')

    # Analytical strains
    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$', '(c)', False)
    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11VertexAnalytical, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$', '(d)', False)
    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$', '(e)', False)
    plot_subfigure(axes[1,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22VertexAnalytical, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$', '(f)', False)
    plot_subfigure(axes[1,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12VertexAnalytical, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$', '(g)', False)
    plot_subfigure(axes[1,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12VertexAnalytical, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$', '(h)', True)

    # MPM
    plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ MPM.', '(i)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain11Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ MPM.', '(j)', False)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$ MPM.', '(k)', False)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain22Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$ MPM.', '(l)', False)
    plot_subfigure(axes[2,4], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$ MPM.', '(m)', False)
    plot_subfigure(axes[2,5], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, strain12Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$ MPM.', '(n)', True)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_map.png",dpi=400)
    #plt.savefig("strain_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()

#---------------------------------------------------------------

if __name__ == "__main__":

    strain_map()
