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

def get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertexx, mpasArray, cmap, vmin, vmax, minX, maxX, minY, maxY):

    patches = []
    colours = []

    minval =  1.0e30
    maxval = -1.0e30

    for iVertex in range(0,nVertices):

        polygonVertices = []

        useVertex = False
        for iCellOnVertex in range(0,vertexDegree[iVertex]):

            iCell = cellsOnVertex[iVertex,iCellOnVertex]

            polygonVertices.append((xCell[iCell],yCell[iCell]))

            if (xCell[iCell] >= minX and xCell[iCell] <= maxX and \
                yCell[iCell] >= minY and yCell[iCell] <= maxY):
                useVertex = True

        if (useVertex and (useVertexx[iVertex] == 1)):
            polygon = Polygon(polygonVertices)
            patches.append(polygon)

            colours.append(mpasArray[iVertex])

            minval = min(minval,mpasArray[iVertex])
            maxval = max(maxval,mpasArray[iVertex])

    patchCollection = PatchCollection(patches, cmap=cmap, rasterized=False)
    patchCollection.set_array(np.array(colours))
    patchCollection.set_linewidth(0)

    patchCollection.set_clim(vmin=vmin,vmax=vmax)

    return patchCollection, minval, maxval

#---------------------------------------------------------------

def plot_subfigure(axes, fig, nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, array, vmin, vmax, minX, maxX, minY, maxY, sciNote=False, diffPlot=False, title=None, subfigureLabel=None, colorbar=True, unityBar=False):

    if (not diffPlot):
        #colourMap = mpl.cm.jet
        colourMap = mpl.cm.seismic
    else:
        #colourMap = mpl.cm.RdBu
        colourMap = mpl.cm.seismic

    patchCollection, minArray, maxArray = get_mpas_patch_collection(nVertices, vertexDegree, cellsOnVertex, xCell, yCell, zCell, useVertex, array, colourMap, vmin, vmax, minX, maxX, minY, maxY)
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

    # grid quad
    fileGrid = Dataset("grid_hex_0082x0094.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    cellsOnCell = fileGrid.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1

    edgesOnCell = fileGrid.variables["edgesOnCell"][:]
    edgesOnCell[:] = edgesOnCell[:] - 1

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]
    latEdge = fileGrid.variables["latEdge"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    xEdge = fileGrid.variables["xEdge"][:]
    yEdge = fileGrid.variables["yEdge"][:]
    zEdge = fileGrid.variables["zEdge"][:]

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    xMax *= 0.2 
    yMax *= 0.2

    fileGrid.close()

    # variational strain cells
    xVar = []
    yVar = []
    zVar = []
    latVar = []

    nEdgesOnCellVar = []
    verticesOnCellVar = []

    useCellVar = []
    
    iVar = 0
    nCellsVar = 0
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

            iVertex = verticesOnCell[iCell,iVertexOnCell]
            iEdgeOnCell1 = iVertexOnCell
            iEdgeOnCell2 = iVertexOnCell - 1
            if (iEdgeOnCell2 < 0): iEdgeOnCell2 = nEdgesOnCell[iCell]-1
            iEdge1 = edgesOnCell[iCell,iEdgeOnCell1]
            iEdge2 = edgesOnCell[iCell,iEdgeOnCell2]

            verticesOnCellVar.append([iVar, iVar+1, iVar+2, iVar+3])

            xVar.append(xCell[iCell])
            yVar.append(yCell[iCell])
            zVar.append(zCell[iCell])
            iVar += 1

            xVar.append(xEdge[iEdge1])
            yVar.append(yEdge[iEdge1])
            zVar.append(zEdge[iEdge1])
            iVar += 1

            xVar.append(xVertex[iVertex])
            yVar.append(yVertex[iVertex])
            zVar.append(zVertex[iVertex])
            iVar += 1

            xVar.append(xEdge[iEdge2])
            yVar.append(yEdge[iEdge2])
            zVar.append(zEdge[iEdge2])
            iVar += 1

            latVar.append(0.25 * (latCell[iCell] + latEdge[iEdge1] + latVertex[iVertex] + latEdge[iEdge2]))

            nEdgesOnCellVar.append(4)

            nCellsVar += 1
            useCellVar.append(1)

    xVar = np.array(xVar)
    yVar = np.array(yVar)
    zVar = np.array(zVar)
    latVar = np.array(latVar)
    nEdgesOnCellVar = np.array(nEdgesOnCellVar)
    verticesOnCellVar = np.array(verticesOnCellVar)
    useCellVar = np.array(useCellVar)

    # ic hex
    fileIC = Dataset("ic_hex_0082x0094.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    strain11VarAnalytical = np.zeros((nCells, maxEdges))
    strain22VarAnalytical = np.zeros((nCells, maxEdges))
    strain12VarAnalytical = np.zeros((nCells, maxEdges))
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            strain11VarAnalytical[iCell,iVertexOnCell] = strain11VertexAnalytical[iVertex]
            strain22VarAnalytical[iCell,iVertexOnCell] = strain22VertexAnalytical[iVertex]
            strain12VarAnalytical[iCell,iVertexOnCell] = strain12VertexAnalytical[iVertex]

    print("Stress divergence: ",
          np.amin(strain11VertexAnalytical), np.amax(strain11VertexAnalytical),
          np.amin(strain22VertexAnalytical), np.amax(strain22VertexAnalytical),
          np.amin(strain12VertexAnalytical), np.amax(strain12VertexAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_hex_wachspress_0082x0094/output.2000.nc","r")

    interiorCell = fileWach.variables["interiorCell"][0,:]

    strain11varWachs = fileWach.variables["strain11var"][0,:]
    strain22varWachs = fileWach.variables["strain22var"][0,:]
    strain12varWachs = fileWach.variables["strain12var"][0,:]

    strain11varWachsDiff = []
    strain22varWachsDiff = []
    strain12varWachsDiff = []
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            strain11varWachsDiff.append(strain11varWachs[iCell,iVertexOnCell] - strain11VarAnalytical[iCell,iVertexOnCell])
            strain22varWachsDiff.append(strain22varWachs[iCell,iVertexOnCell] - strain22VarAnalytical[iCell,iVertexOnCell])
            strain12varWachsDiff.append(strain12varWachs[iCell,iVertexOnCell] - strain12VarAnalytical[iCell,iVertexOnCell])

    print("Wachs: ",
          np.amin(strain11varWachsDiff), np.amax(strain11varWachsDiff),
          np.amin(strain22varWachsDiff), np.amax(strain22varWachsDiff),
          np.amin(strain12varWachsDiff), np.amax(strain12varWachsDiff))

    fileWach.close()

    useVertex = np.ones(nVertices,dtype="i")
    useCell = np.ones(nCells,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            useCell[iCell] = 0
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                useCell[iCell2] = 0
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    # PWL
    filePWL = Dataset("./output_hex_pwl_0082x0094/output.2000.nc","r")

    strain11varPWL = filePWL.variables["strain11var"][0,:]
    strain22varPWL = filePWL.variables["strain22var"][0,:]
    strain12varPWL = filePWL.variables["strain12var"][0,:]

    strain11varPWLDiff = []
    strain22varPWLDiff = []
    strain12varPWLDiff = []
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            strain11varPWLDiff.append(strain11varPWL[iCell,iVertexOnCell] - strain11VarAnalytical[iCell,iVertexOnCell])
            strain22varPWLDiff.append(strain22varPWL[iCell,iVertexOnCell] - strain22VarAnalytical[iCell,iVertexOnCell])
            strain12varPWLDiff.append(strain12varPWL[iCell,iVertexOnCell] - strain12VarAnalytical[iCell,iVertexOnCell])

    print("PWL: ",
          np.amin(strain11varPWLDiff), np.amax(strain11varPWLDiff),
          np.amin(strain22varPWLDiff), np.amax(strain22varPWLDiff),
          np.amin(strain12varPWLDiff), np.amax(strain12varPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./output_hex_weak_0082x0094/output.2000.nc","r")

    strain11weakWeak = fileWeak.variables["strain11weak"][0,:]
    strain22weakWeak = fileWeak.variables["strain22weak"][0,:]
    strain12weakWeak = fileWeak.variables["strain12weak"][0,:]

    strain11weakWeakDiff = (strain11weakWeak - strain11CellAnalytical)
    strain22weakWeakDiff = (strain22weakWeak - strain22CellAnalytical)
    strain12weakWeakDiff = (strain12weakWeak - strain12CellAnalytical)

    print("Weak: ",
          np.amin(strain11weakWeakDiff), np.amax(strain11weakWeakDiff),
          np.amin(strain22weakWeakDiff), np.amax(strain22weakWeakDiff),
          np.amin(strain12weakWeakDiff), np.amax(strain12weakWeakDiff))

    fileWeak.close()


    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    # sinusoid
    #minStrain = -20.0
    #maxStrain =  20.0

    #minStrainDiff = -2.0
    #maxStrainDiff =  2.0

    # linear
    minStrain = None#-1.0
    maxStrain = None# 1.0

    minStrainDiff = None#-1e-2
    maxStrainDiff = None# 1e-2

    fig, axes = plt.subplots(4, 3)

    fig.set_size_inches(7, 6.75)

    # analytical
    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain11VertexAnalytical, minStrain, maxStrain, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ Analytical', r'(a)', True)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain22VertexAnalytical, minStrain, maxStrain, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ Analytical', r'(b)', True)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, useVertex, strain12VertexAnalytical, minStrain, maxStrain, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ Analytical', r'(c)', True)

    # Wachspress
    plot_subfigure(axes[1,0], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, useCellVar, strain11varWachsDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ Wachs.', r'(d)', True)
    plot_subfigure(axes[1,1], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, useCellVar, strain22varWachsDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ Wachs.', r'(e)', True)
    plot_subfigure(axes[1,2], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, useCellVar, strain12varWachsDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ Wachs.', r'(f)', True)

    # PWL
    plot_subfigure(axes[2,0], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, useCellVar, strain11varPWLDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ PWL', r'(g)', True)
    plot_subfigure(axes[2,1], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, useCellVar, strain22varPWLDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ PWL', r'(h)', True)
    plot_subfigure(axes[2,2], fig, nCellsVar, nEdgesOnCellVar, verticesOnCellVar, xVar, yVar, zVar, useCellVar, strain12varPWLDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ PWL', r'(i)', True)

    # Weak
    plot_subfigure(axes[3,0], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, useCell, strain11weakWeakDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{11}$ Weak', r'(j)', True)
    plot_subfigure(axes[3,1], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, useCell, strain22weakWeakDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{22}$ Weak', r'(k)', True)
    plot_subfigure(axes[3,2], fig, nCells, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, useCell, strain12weakWeakDiff, minStrainDiff, maxStrainDiff, xMin, xMax, yMin, yMax, \
                   False, False, r'$\epsilon_{12}$ Weak', r'(l)', True)

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_map.png",dpi=600)
    plt.savefig("strain_map.eps")

    plt.clf()
    plt.cla()
    plt.close()



#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_map()
