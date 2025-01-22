from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

#--------------------------------------------------------

def L2_norm(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

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

def get_norm(filenameIC, filename, latitudeLimit):

    fileIC = Dataset(filenameIC, "r")

    uVelocityAnalytical = fileIC.variables["uVelocity"][:]
    vVelocityAnalytical = fileIC.variables["vVelocity"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    uVelocity = fileMPAS.variables["uVelocity"][0,:]
    vVelocity = fileMPAS.variables["vVelocity"][0,:]

    normU = L2_norm(uVelocity, uVelocityAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    normV = L2_norm(vVelocity, vVelocityAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)

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

def reconstruction_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    tests = ["1-x", "y-z", "latlon", "nonlin"]
    #tests = ["1-x"]

    resolutions = [2562,10242,40962,163842]

    velocities = ['u','v']

    lineColours = ["black","black","gray","silver","blue","blue","green","green","cyan","cyan"]
    lineStyles  = ["-","--","-","--","-","--","-","--","-","--"]

    latitudeLimit = 20.0

    xMin = 6e-3
    xMax = 1e-1

    # quadratic scaling
    scale = 1.0e-3 / math.pow(xMin,2)
    scaleMinQuad = math.pow(xMin,2) * scale
    scaleMaxQuad = math.pow(xMax,2) * scale

    # linear scaling
    scale = 1.0e-3 / math.pow(xMin,1)
    scaleMinLin = math.pow(xMin,1) * scale
    scaleMaxLin = math.pow(xMax,1) * scale


    plt.figure(figsize=(4,3))

    plt.loglog([xMin, xMax],[scaleMinLin,scaleMaxLin],linestyle=':', color='k')
    plt.loglog([xMin, xMax],[scaleMinQuad,scaleMaxQuad],linestyle=':', color='k')

    iPlot = 0

    for test in tests:

        for velocity in velocities:

            x = []
            y = []

            for resolution in resolutions:

                filename = "./output_%s_%i/output.2000.nc" %(test,resolution)
                filenameIC = "./ic_%s_%i.nc" %(test,resolution)

                print(filename, filenameIC)

                normU, normV = get_norm(filenameIC, filename, latitudeLimit)

                x.append(get_resolution(filename, latitudeLimit))

                if (velocity == "u"):
                    y.append(normU)
                elif (velocity == "v"):
                    y.append(normV)

            plt.loglog(x,y, marker='o', color=lineColours[iPlot], ls=lineStyles[iPlot], markersize=5.0)
            print(velocity, x, y)

            iPlot = iPlot + 1

    legendLabels = ["linear scaling", "quadratic scaling", "1", "x", "y", "z", "lat", "lon", "quadratic", "1/cos(lat)+tan(lat)"]
    plt.legend(legendLabels, frameon=False, loc=3, fontsize=6, handlelength=4)

    ax = plt.gca()
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel(r"$L_2$ error norm")
    ax.set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("reconstruction_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    reconstruction_scaling()
