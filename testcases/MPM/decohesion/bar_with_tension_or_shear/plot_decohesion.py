from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
import glob
import os

#-------------------------------------------------------------------------------

def plot_decohesion(gridFilename):

    if (not os.path.isdir("plots")):
                 os.mkdir("plots")

    filein = Dataset(gridFilename,"r")
    nCells = len(filein.dimensions["nCells"])
    nEdges = len(filein.dimensions["nEdges"])
    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    verticesOnCell = filein.variables["verticesOnCell"][:]-1
    edgesOnCell = filein.variables["edgesOnCell"][:]-1
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

    filenames = sorted(glob.glob("./output/particles_output.*"))
    nTimes = len(filenames)

    for iTime in range(0, nTimes):
    #iTime = 2
    #if (iTime == 2):

        filein = Dataset(filenames[iTime],"r")
        nParticles = len(filein.dimensions["nParticles"])
        posnMP = filein.variables["posnMP"][:,:,:]
        decohesionOpeningMP = filein.variables["decohesionOpeningMP"][:,:,:]
        decohesionCrackEndPoints = filein.variables["decohesionCrackEndPoints"][:,:,:]
        decohesionCrackEdgeIndex = filein.variables["decohesionCrackEdgeIndex"][:,:,:]
        filein.close()

        # plot decohesionOpeningMP at posnMP
        fig, axis = plt.subplots(figsize=(10,10))
        color = abs(decohesionOpeningMP[0,:,0]) + abs(decohesionOpeningMP[0,:,1])
        scatter = axis.scatter(posnMP[0,:,0], posnMP[0,:,1], c = color)
        axis.set_xlim(xmin,xmax)
        axis.set_ylim(ymin,ymax)
        axis.set_xlabel("x")
        axis.set_ylabel("y")

        plt.colorbar(scatter)
        plt.tight_layout()
        filenameOut = "plots/decohesion_%i.png"%(iTime)
        plt.savefig(filenameOut,dpi=600)
        plt.close()

        # plot the mesh and cracks
        patchesCell = []
        for iCell in range(0,nCells):
           vertices = []
           for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
               iVertex = verticesOnCell[iCell,iVertexOnCell]
               vertices.append((xVertex[iVertex],yVertex[iVertex]))
           patchesCell.append(Polygon(vertices,closed=True, edgecolor="teal", fill=False, linewidth=0.1))

        pc = PatchCollection(patchesCell, match_original=True)

        crackLines = []
        for iCell in range(0, nCells):
           crackIndex = decohesionCrackEdgeIndex[0,iCell,:] - 1
           if (crackIndex[0] >= 0 and crackIndex[1] >=0):
              iEdges = edgesOnCell[iCell, crackIndex]
              crackLines.append([(decohesionCrackEndPoints[0,iEdges[0],0],decohesionCrackEndPoints[0,iEdges[0],1]),
                                 (decohesionCrackEndPoints[0,iEdges[1],0],decohesionCrackEndPoints[0,iEdges[1],1])])

        lc = LineCollection(crackLines)

        fig, axis = plt.subplots(figsize=(10,10))
        axis.add_collection(pc)
        axis.add_collection(lc)
        axis.scatter(decohesionCrackEndPoints[0,:,0], decohesionCrackEndPoints[0,:,1], marker='o')
        axis.set_xlim((xmin-0.05*lxy,xmax+0.05*lxy))
        axis.set_ylim((ymin-0.05*lxy,ymax+0.05*lxy))
        axis.set_aspect('equal')
        axis.set_xlabel("x")
        axis.set_ylabel("y")

        plt.tight_layout()
        filenameOut = "plots/cracks_%i.png"%(iTime)
        plt.savefig(filenameOut,dpi=600)
        plt.close()

#-------------------------------------------------------------------------------
