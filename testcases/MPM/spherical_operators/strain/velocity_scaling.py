from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import numpy as np

#--------------------------------------------------------

def L2_norm_vertex(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (math.fabs(latVertex[iVertex]) > latitudeLimit * degreesToRadians):

            norm  = norm  + areaTriangle[iVertex] * math.pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * math.pow(analytical[iVertex],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_vertex(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    uVelocityAnalytical = fileIC.variables["uVelocity"][:]
    vVelocityAnalytical = fileIC.variables["vVelocity"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])
    nCells = len(fileMPAS.dimensions["nCells"])
    vertexDegree = len(fileMPAS.dimensions["vertexDegree"])

    latVertex = fileMPAS.variables["latVertex"][:]
    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    uVelocity= fileMPAS.variables["uVelocity"][:]
    vVelocity= fileMPAS.variables["vVelocity"][:]

    fileMPAS.close()

    normU = L2_norm_vertex(uVelocity[0,:], uVelocityAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normV = L2_norm_vertex(vVelocity[0,:], vVelocityAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

    return normU, normV

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

    uvVelAnalyticalMP = fileIC.variables["uvVelMP"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    latParticle = fileMPAS.variables["latCellMP"][:]

    uvVelMP = fileMPAS.variables["uvVelMP"][:]

    normU = L2_norm_particle(uvVelMP[0,:,0], uvVelAnalyticalMP[:,0], nParticles, latParticle[0,:], latitudeLimit)
    normV = L2_norm_particle(uvVelMP[0,:,1], uvVelAnalyticalMP[:,1], nParticles, latParticle[0,:], latitudeLimit)

    fileMPAS.close()

    return normU, normV

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

def velocity_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    velocities = ["u","v"]

    resolutions = [2562,10242,40962,163842]

    methods = ["mpmvar"]

    lineColours = ["black","grey","red"]

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    latitudeLimit = 20.0

    fig, axes = plt.subplots(2,1,figsize=(5,10))

    iVelo = 0
    for velocity in velocities:

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

        axes[iVelo].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes[iVelo].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

        iPlot = 0
        for method in methods:

            x = []
            y = []

            for resolution in resolutions:

                filenames = sorted(glob.glob("./output_%s_%i/particles_output*" %(method,resolution)))
                filename = filenames[-1]
                filenameIC = "particles_%s.nc" %resolution

                normU, normV = get_norm_particle(filenameIC, filename, latitudeLimit)

                filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                x.append(get_resolution(filename, latitudeLimit))

                if (velocity == "u"):
                    y.append(normU)
                elif (velocity == "v"):
                    y.append(normV)

            iPlot = iPlot + 1

        axes[iVelo].loglog(x,y, marker='o', color=lineColours[1], ls="solid", markersize=5.0, label="VelMP" )
        print('resolution',x)
        print('normMP',y)

        iPlot = 0
        for method in methods:

            x = []
            y = []

            for resolution in resolutions:

                filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                filenameIC = "./ic_%i.nc" %(resolution)

                normU, normV = get_norm_vertex(filenameIC, filename, latitudeLimit)

                x.append(get_resolution(filename, latitudeLimit))

                if (velocity == "u"):
                    y.append(normU)
                elif (velocity == "v"):
                    y.append(normV)

            iPlot = iPlot + 1

        axes[iVelo].loglog(x,y, marker='o', color=lineColours[2], ls="solid", markersize=5.0, label="VelVertex" )
        print('resolution',x)
        print('normVertex',y)


        axes[iVelo].legend(frameon=False, loc=2, fontsize=8, handlelength=4)

        axes[iVelo].set_xlabel("Grid resolution")
        axes[iVelo].set_ylabel(r"$L_2$ error norm")

        iVelo = iVelo + 1


    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("velocity_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    velocity_scaling()
