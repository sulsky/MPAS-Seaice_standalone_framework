from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import glob

#-------------------------------------------------------------------------------

def plot_particles_scatter():

    arrayName = "iceVolumeCellMP"

    ymin = 200000.0
    ymax = 300000.0



    filenames = sorted(glob.glob("./output_particles/particles_output.*"))

    nTimes = len(filenames)

    fig, axis = plt.subplots()

    cmap = plt.cm.jet
    colors = cmap(np.linspace(0, 1, nTimes))

    for iTime in range(0,nTimes):

        filein = Dataset(filenames[iTime],"r")

        try:
            nParticles = len(filein.dimensions["nParticles"])
        except:
            continue

        if (nParticles > 0):

            #print(filenames[iTime])

            x = filein.variables["posnMP"][0,:,0]
            y = filein.variables["posnMP"][0,:,1]

            indexToCellIDMP = filein.variables["indexToCellIDMP"][0,:]

            indices = np.where((y >= ymin) & (y <= ymax))[0]

            arrayParticle = filein.variables[arrayName][0,:]

            for iParticle in range(0,nParticles):
                if (arrayParticle[iParticle] > 1.0):
                    print(iParticle,arrayParticle[iParticle],indexToCellIDMP[iParticle])

            filein.close()

            axis.scatter(x[indices], arrayParticle[indices], color=colors[iTime], s=1, rasterized=True)

    axis.set_xlim((0,1000000))
    axis.set_ylim((0,0.6))
    plt.tight_layout()
    plt.savefig("particles_scatter.png",dpi=300)
    plt.savefig("particles_scatter.pdf",dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_particles_scatter()
