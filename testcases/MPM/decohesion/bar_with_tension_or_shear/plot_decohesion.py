from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
import argparse
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

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='gridFilename')
    args = parser.parse_args()
    plot_decohesion(args.gridFilename)
