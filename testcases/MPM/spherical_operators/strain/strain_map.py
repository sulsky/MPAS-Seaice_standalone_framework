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

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if (colorbar):
        cb = fig.colorbar(patchCollection,cax=cax)
        if (unityBar):
            cb.ax.set_yticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
        if (sciNote):
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()
    else:
        cax.set_axis_off()

#---------------------------------------------------------------

def scatter_plot(axes, fig, posnMP, color, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    scatter = axes.scatter(posnMP[:,0], posnMP[:,1], posnMP[:,2], c = color, cmap=colourMap)
    axes.set_ylim([minY,maxY])
    axes.set_xlim([minX,maxX])
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_aspect('equal', adjustable='box')

    if (title != None):
        axes.set_title(title, fontsize=8)

    if (subfigureLabel != None):
        axes.text(0.02, 0.89, subfigureLabel, verticalalignment='bottom', horizontalalignment='left',transform=axes.transAxes, fontsize=8)

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    scatter.set_clim(vmin=vmin, vmax=vmax)
    if (colorbar):
        cb = fig.colorbar(scatter,cax=cax)
        if (unityBar):
            cb.ax.set_yticklabels(['-1.0','-0.5','0.0','0.5','1.0'])
        if (sciNote):
            cb.formatter.set_powerlimits((0, 0))
            cb.update_ticks()
    else:
        cax.set_axis_off()

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

    fileIC = Dataset("particles_40962.nc","r")

    strain11AnalyticalMP = fileIC.variables["strainAnalyticalMP"][:,0]
    strain22AnalyticalMP = fileIC.variables["strainAnalyticalMP"][:,1]
    strain12AnalyticalMP = fileIC.variables["strainAnalyticalMP"][:,2]
    posnMP = fileIC.variables["posnMP"][:,:]

    fileIC.close()

    fileMPM = Dataset("./output_mpm_40962/particles_output.2000-01-01_01.00.00.nc","r")

    strain11 = fileMPM.variables["strainRateMP"][0,:,0]
    strain22 = fileMPM.variables["strainRateMP"][0,:,1]
    strain12 = fileMPM.variables["strainRateMP"][0,:,2]

    strain11Diff = strain11 - strain11AnalyticalMP
    strain22Diff = strain22 - strain22AnalyticalMP
    strain12Diff = strain12 - strain12AnalyticalMP

    print("MPM: ",
          np.amin(strain11Diff), np.amax(strain11Diff),
          np.amin(strain22Diff), np.amax(strain22Diff),
          np.amin(strain12Diff), np.amax(strain12Diff))

    fileMPM.close()

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(4, 6)
    fig.set_size_inches(9, 6)

    minStrain = -3.5
    maxStrain =  3.5

    minDiff = -0.10
    maxDiff =  0.10

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
    scatter_plot(axes[2,0], fig, posnMP, strain11, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{11}$ MPM.', '(i)', False)
    scatter_plot(axes[2,1], fig, posnMP, strain11, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{11}$ MPM.', '(j)', False)
    scatter_plot(axes[2,2], fig, posnMP, strain22, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{22}$ MPM.', '(k)', False)
    scatter_plot(axes[2,3], fig, posnMP, strain22, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{22}$ MPM.', '(l)', False)
    scatter_plot(axes[2,4], fig, posnMP, strain12, minStrain, maxStrain, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$\epsilon_{12}$ MPM.', '(m)', False)
    scatter_plot(axes[2,5], fig, posnMP, strain12, minStrain, maxStrain, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'$\epsilon_{12}$ MPM.', '(n)', True)

    scatter_plot(axes[3,0], fig, posnMP, strain11Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\epsilon_{11}$ MPM.', '(i)', False)
    scatter_plot(axes[3,1], fig, posnMP, strain11Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\epsilon_{11}$ MPM.', '(j)', False)
    scatter_plot(axes[3,2], fig, posnMP, strain22Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\epsilon_{22}$ MPM.', '(k)', False)
    scatter_plot(axes[3,3], fig, posnMP, strain22Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\epsilon_{22}$ MPM.', '(l)', False)
    scatter_plot(axes[3,4], fig, posnMP, strain12Diff, minDiff, maxDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, True, r'$\epsilon_{12}$ MPM.', '(m)', False)
    scatter_plot(axes[3,5], fig, posnMP, strain12Diff, minDiff, maxDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, True, r'$\epsilon_{12}$ MPM.', '(n)', True)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_map.png",dpi=400)
    #plt.savefig("strain_map_3.png", bbox_inches="tight",dpi=2000)

    plt.clf()
    plt.cla()
    plt.close()

#---------------------------------------------------------------

if __name__ == "__main__":

    strain_map()
