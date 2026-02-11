from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
import glob
import os

#-------------------------------------------------------------------------------

def plot_disks():

    if (not os.path.isdir("plots")):
                 os.mkdir("plots")

    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0

    filenames = sorted(glob.glob("./output/particles_output.*"))
    nTimes = len(filenames)

    for iTime in range(0, nTimes):
    #iTime = 2
    #if (iTime == 2):

        filein = Dataset(filenames[iTime],"r")
        nParticles = len(filein.dimensions["nParticles"])
        posnMP = filein.variables["posnMP"][:,:,:]
        uvVelMP = filein.variables["uvVelMP"][:,:,:]
        filein.close()

        # plot decohesionOpeningMP at posnMP
        fig, axis = plt.subplots(figsize=(10,10))
        scatter = axis.scatter(posnMP[0,:,0], posnMP[0,:,1], c = uvVelMP[0,:,1])
        axis.set_xlim(xmin,xmax)
        axis.set_ylim(ymin,ymax)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_aspect('equal')

        plt.colorbar(scatter)
        plt.tight_layout()
        filenameOut = "plots/disks_%i.png"%(iTime)
        plt.savefig(filenameOut,dpi=600)
        plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_disks()

