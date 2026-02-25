from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import numpy as np

#--------------------------------------------------------

def L2_norm_particle(numerical, analytical, nParticles, lat, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iParticle in range(0,nParticles):

        if (math.fabs(lat[iParticle]) > latitudeLimit * degreesToRadians):

            norm  = norm + math.pow(numerical[iParticle] - analytical[iParticle],2)

            denom = denom + math.pow(analytical[iParticle],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_particle(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    strainAnalyticalMP = fileIC.variables["strainAnalyticalMP"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    latParticle = fileMPAS.variables["latCellMP"][:]

    strainRateMP = fileMPAS.variables["strainRateMP"][:]

    normE11 = L2_norm_particle(strainRateMP[0,:,0], strainAnalyticalMP[:,0], nParticles, latParticle[0,:], latitudeLimit)
    normE22 = L2_norm_particle(strainRateMP[0,:,1], strainAnalyticalMP[:,1], nParticles, latParticle[0,:], latitudeLimit)
    normE12 = L2_norm_particle(strainRateMP[0,:,2], strainAnalyticalMP[:,2], nParticles, latParticle[0,:], latitudeLimit)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def get_resolution(filename, latitudeLimit):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    degreesToRadians = math.pi / 180.0

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        if (math.fabs(latEdge[iEdge]) > latitudeLimit * degreesToRadians):
            resolution = resolution + dcEdge[iEdge]
            denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

#--------------------------------------------------------

def scaling_lines(axis, xMin, xMax, yMin):

    # linear scaling
    scale = yMin / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale

    # quadratic scaling
    scale = yMin / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    axis.loglog([xMin, xMax], [scaleMinLin,  scaleMaxLin],  linestyle=':', color='k')
    axis.loglog([xMin, xMax], [scaleMinQuad, scaleMaxQuad], linestyle=':', color='k')

#--------------------------------------------------------

def strain_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    strains = ["strain11","strain22","strain12"]

    resolutions = [2562,10242,40962,163842]

    methods = ["mpm"]

    lineColours = ["black","grey","red"]

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    labels = [r"$\varepsilon_{11}$", r"$\varepsilon_{22}$", r"$\varepsilon_{12}$"]

    latitudeLimit = 20.0

    fig, axes = plt.subplots(figsize=(3,3))

    iPlot = 0
    for strain in strains:

        xMin = 4e-2
        xMax = 8e-2

        # linear scaling
        scale = 1e-3 / math.pow(xMin,1)
        scaleMinLin = math.pow(xMin,1) * scale
        scaleMaxLin = math.pow(xMax,1) * scale

        # quadratic scaling
        scale = 1e-3 / math.pow(xMin,2)
        scaleMinQuad = math.pow(xMin,2) * scale
        scaleMaxQuad = math.pow(xMax,2) * scale

        axes.loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes.loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')


        for method in methods:

            x = []
            y = []

            for resolution in resolutions:

                filenames = sorted(glob.glob("./output_%s_%i/particles_output*" %(method,resolution)))
                filename = filenames[-1]
                filenameIC = "particles_%s.nc" %resolution

                normE11, normE22, normE12 = get_norm_particle(filenameIC, filename, latitudeLimit)

                filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                x.append(get_resolution(filename, latitudeLimit))
                if (strain == "strain11"):
                    y.append(normE11)
                elif (strain == "strain22"):
                    y.append(normE22)
                elif (strain == "strain12"):
                    y.append(normE12)

            axes.loglog(x,y, marker='o', color=lineColours[iPlot], ls="solid", markersize=5.0, label=labels[iPlot])

        iPlot = iPlot + 1

        axes.legend(frameon=False, loc=2, fontsize=8, handlelength=4)

        axes.set_xlabel("Grid resolution")
        axes.set_ylabel(r"$L_2$ error norm")
        #axes[iStrain].set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("strain_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_scaling()
