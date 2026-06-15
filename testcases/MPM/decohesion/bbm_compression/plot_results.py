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

def plot_results(gridFilename):

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
    iTime = -1
    filein = Dataset(filenames[iTime],"r")
    nParticles = len(filein.dimensions["nParticles"])
    posnMP = filein.variables["posnMP"][:,:,:]
    uvVelMP = filein.variables["uvVelMP"][:,:,:]
    decohesionOpeningMP = filein.variables["decohesionOpeningMP"][:,:,:]
    filein.close()

    # plot decohesionOpeningMP at posnMP
    fig, axes = plt.subplots(2,2)
    color = decohesionOpeningMP[0,:,0]
    axes[0,0].set_xlim(xmin,xmax)
    axes[0,0].set_ylim(ymin,ymax)
    M = axes[0,0].transData.get_matrix()
    xscale = M[0,0]
    yscale = M[1,1]
    size = xscale * yscale * lx * ly / nCells / 36
    scatter = axes[0,0].scatter(posnMP[0,:,0], posnMP[0,:,1], c = color, s = size)
    axes[0,0].set_xlabel("x")
    axes[0,0].set_ylabel("y")
    axes[0,0].set_title("normal opening")
    scatter.set_clim(vmin=0, vmax=10)
    fig.colorbar(scatter)
    axes[0,0].set_aspect('equal')

    color = decohesionOpeningMP[0,:,1]
    axes[0,1].set_xlim(xmin,xmax)
    axes[0,1].set_ylim(ymin,ymax)
    M = axes[0,1].transData.get_matrix()
    xscale = M[0,0]
    yscale = M[1,1]
    size = xscale * yscale * lx * ly / nCells / 36
    scatter = axes[0,1].scatter(posnMP[0,:,0], posnMP[0,:,1], c = color, s = size)
    axes[0,1].set_xlabel("x")
    axes[0,1].set_ylabel("y")
    axes[0,1].set_title("tangential opening")
    scatter.set_clim(vmin=-10, vmax=10)
    fig.colorbar(scatter)
    axes[0,1].set_aspect('equal')

    color = abs(decohesionOpeningMP[0,:,0]) + abs(decohesionOpeningMP[0,:,1])
    axes[1,0].set_xlim(xmin,xmax)
    axes[1,0].set_ylim(ymin,ymax)
    M = axes[1,0].transData.get_matrix()
    xscale = M[0,0]
    yscale = M[1,1]
    size = xscale * yscale * lx * ly / nCells / 36
    scatter = axes[1,0].scatter(posnMP[0,:,0], posnMP[0,:,1], c = color, s = size)
    axes[1,0].set_xlabel("x")
    axes[1,0].set_ylabel("y")
    axes[1,0].set_title("combined (1-norm) opening")
    scatter.set_clim(vmin=0, vmax=10)
    fig.colorbar(scatter)
    axes[1,0].set_aspect('equal')

    color = uvVelMP[0,:,0]
    axes[1,1].set_xlim(xmin,xmax)
    axes[1,1].set_ylim(ymin,ymax)
    M = axes[1,1].transData.get_matrix()
    xscale = M[0,0]
    yscale = M[1,1]
    size = xscale * yscale * lx * ly / nCells / 36
    scatter = axes[1,1].scatter(posnMP[0,:,0], posnMP[0,:,1], c = color, s = size)
    axes[1,1].set_xlabel("x")
    axes[1,1].set_ylabel("y")
    axes[1,1].set_title("u-Velocity")
    scatter.set_clim(vmin=-0.05, vmax=1.e-5)
    fig.colorbar(scatter)
    axes[1,1].set_aspect('equal')

    plt.tight_layout()
    filenameOut = "decohesion.png"
    plt.savefig(filenameOut,dpi=600)
    plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='gridFilename')
    args = parser.parse_args()
    plot_results(args.gridFilename)
