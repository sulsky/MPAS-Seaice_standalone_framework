from netCDF4 import Dataset
import math
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

#--------------------------------------------------------

def L2_norm(numerical, analytical, nPoints, area, lat, latitudeLimit):

    norm  = 0.0
    denom = 0.0

    degreesToRadians = math.pi / 180.0

    for iPoint in range(0,nPoints):

        if (math.fabs(lat[iPoint]) > latitudeLimit * degreesToRadians):

            norm  += area[iPoint] * math.pow(numerical[iPoint] - analytical[iPoint],2)

            denom += area[iPoint] * math.pow(analytical[iPoint],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_area(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    analytical = fileIC.variables["iceAreaCellMP"][:]

    areaMP = fileIC.variables["areaMP"][:]

    fileIC.close()


    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    latCellMP = fileMPAS.variables["latCellMP"][0,:]

    numerical = fileMPAS.variables["windSpeedMP"][0,:]

    fileMPAS.close()


    norm = L2_norm(numerical, analytical, nParticles, areaMP, latCellMP, latitudeLimit)

    return norm

#--------------------------------------------------------

def get_resolution(filenameMesh, latitudeLimit):

    fileMPAS = Dataset(filenameMesh, "r")

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

def interpolation_cell_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    tests = ["1", "1+x", "1+y", "1+z", "lat", "lon", "nonlin"]

    resolutions = [2562,10242,40962,163842]

    latitudeLimit = 0.0

    xMin = 6e-3
    xMax = 1e-1

    # quadratic scaling
    scale = 1.0e-4 / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    # linear scaling
    scale = 1.0e-4 / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale


    plt.figure(figsize=(4,3))

    plt.loglog([xMin, xMax],[scaleMinLin,scaleMaxLin],linestyle='--', color='k')
    plt.loglog([xMin, xMax],[scaleMinQuad,scaleMaxQuad],linestyle=':', color='k')

    lineColours = ["black","red","blue","green","cyan","gray","green"]
    lineStyles  = ["-","-","--",":","-","-","-"]

    iPlot = 0

    for test in tests:

        x = []
        y = []

        for resolution in resolutions:

            outputDir = "./output_%s_%i/" %(test,resolution)
            if (not os.path.isdir(outputDir)):
                raise Exception("Missing output directory")
            filenames = sorted(glob.glob(outputDir+"/particles_output*"))
            if (len(filenames) == 0):
                raise Exception("No output particle files: "+outputDir+"/particles_output*")

            filename = filenames[0]
            filenameIC = "./particles_%s_%i.nc" %(test,resolution)

            print(filename, filenameIC)
            if (not os.path.exists(filename)):
                raise Exception("Missing output file: %s" %(filename))
            if (not os.path.exists(filenameIC)):
                raise Exception("Missing IC file: %s" %(filenameIC))

            filenameMesh = "./grid.%i.nc" %(resolution)
            x.append(get_resolution(filenameMesh, latitudeLimit))

            norm = get_norm_area(filenameIC, filename, latitudeLimit)
            y.append(norm)

        plt.loglog(x,y, marker='o', color=lineColours[iPlot], ls=lineStyles[iPlot], markersize=5.0)
        print(x, y)

        iPlot = iPlot + 1

    legendLabels = ["linear scaling", "quadratic scaling", "1", "1+x", "1+y", "1+z", "lat", "lon", "2+cos(lat)cos(lat)cos(2lon)"]
    plt.legend(legendLabels, frameon=False, loc=4, fontsize=6, handlelength=4)

    ax = plt.gca()
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel(r"$L_2$ error norm")
    ax.set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("interpolation_cell_scaling.png",dpi=400)
    plt.savefig("interpolation_cell_scaling.pdf")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    interpolation_cell_scaling()
