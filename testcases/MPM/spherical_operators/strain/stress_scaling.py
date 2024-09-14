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

    stress11VertexAnalytical = fileIC.variables["stress11VertexAnalytical"][:]
    stress22VertexAnalytical = fileIC.variables["stress22VertexAnalytical"][:]
    stress12VertexAnalytical = fileIC.variables["stress12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])
    nCells = len(fileMPAS.dimensions["nCells"])
    vertexDegree = len(fileMPAS.dimensions["vertexDegree"])

    latVertex = fileMPAS.variables["latVertex"][:]
    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    areaCell = fileMPAS.variables["areaCell"][:]
    cellsOnVertex = fileMPAS.variables["cellsOnVertex"][:]
    cellVerticesAtVertex = fileMPAS.variables["cellVerticesAtVertex"][:]

    stress11var = fileMPAS.variables["stress11var"][0,:]
    stress22var = fileMPAS.variables["stress22var"][0,:]
    stress12var = fileMPAS.variables["stress12var"][0,:]

    fileMPAS.close()

    stress11varAvg, stress22varAvg, stress12varAvg = average_stress(stress11var, stress22var, stress12var, nVertices, nCells, vertexDegree, cellsOnVertex, cellVerticesAtVertex, areaCell)

    normE11 = L2_norm_vertex(stress11varAvg, stress11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normE22 = L2_norm_vertex(stress22varAvg, stress22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normE12 = L2_norm_vertex(stress12varAvg, stress12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

    return normE11, normE22, normE12

#--------------------------------------------------------

def average_stress(stress11var, stress22var, stress12var, nVertices, nCells, vertexDegree, cellsOnVertex, cellVerticesAtVertex, areaCell):

    stress11VarAvg = np.zeros([nVertices])
    stress22VarAvg = np.zeros([nVertices])
    stress12VarAvg = np.zeros([nVertices])

    for iVertex in range(0, nVertices):

          stress11avg = 0.0
          stress22avg = 0.0
          stress12avg = 0.0
          denominator = 0.0

          for iVertexDegree in range(0, vertexDegree):

             iCell = cellsOnVertex[iVertex, iVertexDegree] -1

             if (iCell < nCells):

                iVertexOnCell = cellVerticesAtVertex[iVertex, iVertexDegree]

                stress11avg = stress11avg + stress11var[iCell, iVertexOnCell - 1] * areaCell[iCell]
                stress22avg = stress22avg + stress22var[iCell, iVertexOnCell - 1] * areaCell[iCell]
                stress12avg = stress12avg + stress12var[iCell, iVertexOnCell - 1] * areaCell[iCell]
                denominator = denominator + areaCell[iCell]


          stress11VarAvg[iVertex] = stress11avg / denominator
          stress22VarAvg[iVertex] = stress22avg / denominator
          stress12VarAvg[iVertex] = stress12avg / denominator

    return stress11VarAvg, stress22VarAvg, stress12VarAvg

#--------------------------------------------------------

def L2_norm_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    degreesToRadians = math.pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (math.fabs(latCell[iCell]) > latitudeLimit * degreesToRadians):

            norm  = norm  + areaCell[iCell] * math.pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * math.pow(analytical[iCell],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_cell(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    stress11CellAnalytical = fileIC.variables["stress11CellAnalytical"][:]
    stress22CellAnalytical = fileIC.variables["stress22CellAnalytical"][:]
    stress12CellAnalytical = fileIC.variables["stress12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    latCell = fileMPAS.variables["latCell"][:]

    areaCell = fileMPAS.variables["areaCell"][:]

    stress11weak = fileMPAS.variables["stress11weak"][0,:]
    stress22weak = fileMPAS.variables["stress22weak"][0,:]
    stress12weak = fileMPAS.variables["stress12weak"][0,:]

    normE11 = L2_norm_cell(stress11weak, stress11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    normE22 = L2_norm_cell(stress22weak, stress22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    normE12 = L2_norm_cell(stress12weak, stress12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    fileMPAS.close()

    return normE11, normE22, normE12

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

    stressAnalyticalMP = fileIC.variables["stressAnalyticalMP"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    latLonParticle = fileMPAS.variables["posnLatLonGeoMP"][:]

    stressMP = fileMPAS.variables["stressMP"][:]

    normE11 = L2_norm_particle(stressMP[0,:,0], stressAnalyticalMP[:,0], nParticles, latLonParticle[0,:,0], latitudeLimit)
    normE22 = L2_norm_particle(stressMP[0,:,1], stressAnalyticalMP[:,1], nParticles, latLonParticle[0,:,0], latitudeLimit)
    normE12 = L2_norm_particle(stressMP[0,:,2], stressAnalyticalMP[:,2], nParticles, latLonParticle[0,:,0], latitudeLimit)

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

def stress_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    stresses = ["stress11","stress22","stress12"]

    resolutions = [2562,10242,40962,163842]

    methods = ["mpmvar","mpmweak"]

    lineColours = ["black","grey","red","blue","green"]

    lineStyles = {"hex":"solid",
                  "quad":"dashed"}

    latitudeLimit = 20.0

    fig, axes = plt.subplots(3,1,figsize=(5,10))

    iStress = 0
    for stress in stresses:

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

        axes[iStress].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes[iStress].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

        iPlot = 0
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

                if (stress == "stress11"):
                    y.append(normE11)
                elif (stress == "stress22"):
                    y.append(normE22)
                elif (stress == "stress12"):
                    y.append(normE12)


            axes[iStress].loglog(x,y, marker='o', color=lineColours[iPlot], ls="solid", markersize=5.0, label="StressMP_%s" %(method) )

            iPlot = iPlot + 1

        axes[iStress].legend(frameon=False, loc=2, fontsize=8, handlelength=4)

        axes[iStress].set_xlabel("Grid resolution")
        axes[iStress].set_ylabel(r"$L_2$ error norm")

        iStress = iStress + 1

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_scalingMP.png", dpi=400)


    fig, axes = plt.subplots(3,1,figsize=(5,10))

    iStress = 0
    for stress in stresses:

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

        axes[iStress].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
        axes[iStress].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

        iPlot = 0
        for method in methods:

            x = []
            y = []

            for resolution in resolutions:

                filename = "./output_%s_%i/output.2000.nc" %(method,resolution)
                filenameIC = "./ic_%i.nc" %(resolution)

                print(filename, filenameIC)

                if (method == "mpmvar"):
                   normE11, normE22, normE12 = get_norm_vertex(filenameIC, filename, latitudeLimit)
                else:
                   normE11, normE22, normE12 = get_norm_cell(filenameIC, filename, latitudeLimit)

                x.append(get_resolution(filename, latitudeLimit))

                if (stress == "stress11"):
                    y.append(normE11)
                elif (stress == "stress22"):
                    y.append(normE22)
                elif (stress == "stress12"):
                    y.append(normE12)


            axes[iStress].loglog(x,y, marker='o', color=lineColours[iPlot], ls="solid", markersize=5.0, label="stress_%s" %(method))

            iPlot = iPlot + 1


        axes[iStress].legend(frameon=False, loc=2, fontsize=8, handlelength=4)

        axes[iStress].set_xlabel("Grid resolution")
        axes[iStress].set_ylabel(r"$L_2$ error norm")

        iStress = iStress + 1


    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("stress_scaling.png", dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_scaling()
