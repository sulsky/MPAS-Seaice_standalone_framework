from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import os

#-------------------------------------------------------------------------------

def plot_mesh(gridFilename):

    filein = Dataset(gridFilename,"r")
    nCells = len(filein.dimensions["nCells"])
    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    verticesOnCell = filein.variables["verticesOnCell"][:]-1
    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]
    filein.close()

    xmin = np.amin(xVertex)
    xmax = np.amax(xVertex)
    ymin = np.amin(yVertex)
    ymax = np.amax(yVertex)
    lx = xmax - xmin
    ly = ymax - ymin
    lxy = max(lx,ly)

    patchesCell = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex],yVertex[iVertex]))
        patchesCell.append(Polygon(vertices,closed=True, edgecolor="teal", fill=False, linewidth=0.1))

    pc = PatchCollection(patchesCell, match_original=True)

    fig, axes = plt.subplots()

    axes.add_collection(pc)
    axes.set_xlim((xmin-0.05*lxy,xmax+0.05*lxy))
    axes.set_ylim((ymin-0.05*lxy,ymax+0.05*lxy))
    axes.set_aspect('equal')

    plt.tight_layout()
    filenameOut = os.path.splitext(os.path.basename(gridFilename))[0]+".png"
    plt.savefig(filenameOut,dpi=600)

#-------------------------------------------------------------------------------
