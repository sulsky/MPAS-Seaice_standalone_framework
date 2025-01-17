from netCDF4 import Dataset
import math
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl

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

    normU = L2_norm_particle(uvVelMP[0,:,0], 2.0*uvVelAnalyticalMP[:,0], nParticles, latParticle[0,:], latitudeLimit)
    normV = L2_norm_particle(uvVelMP[0,:,1], 2.0*uvVelAnalyticalMP[:,1], nParticles, latParticle[0,:], latitudeLimit)

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

def interpolation_scaling():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    tests = ["1-x", "y-z", "latlon", "nonlin"]
    #tests = ["1-x"]

    resolutions = [2562,10242,40962,163842]

    velocities = ['u','v']

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

    lineColours = ["black","black","red","red","blue","blue","green","green","cyan","cyan"]
    lineStyles  = ["-","--","-","--","-","--","-","--","-","--"]

    iPlot = 0

    for test in tests:

        for velocity in velocities:

            x = []
            y = []

            for resolution in resolutions:

                filenames = sorted(glob.glob("./output_%s_%i/particles_output*" %(test,resolution)))
                filename = filenames[-1]
                filenameIC = "particles_%s_%s.nc" %(test,resolution)

                print(filename, filenameIC)

                normU, normV = get_norm_particle(filenameIC, filename, latitudeLimit)

                filename = "./output_%s_%i/output.2000.nc" %(test,resolution)
                x.append(get_resolution(filename, latitudeLimit))

                if (velocity == "u"):
                    y.append(normU)
                elif (velocity == "v"):
                    y.append(normV)

            plt.loglog(x,y, marker='o', color=lineColours[iPlot], ls=lineStyles[iPlot], markersize=5.0)
            print(velocity, x, y)

            iPlot = iPlot + 1


    legendLabels = ["linear scaling", "quadratic scaling", "1", "x", "y", "z", "lat", "lon", "quadratic", "1/cos(lat)+tan(lat)"]
    #legendLabels = ["linear scaling", "quadratic scaling", "1", "x"]
    plt.legend(legendLabels, frameon=False, loc=4, fontsize=6, handlelength=4)

    ax = plt.gca()
    ax.set_xlabel("Grid resolution")
    ax.set_ylabel(r"$L_2$ error norm")
    ax.set_xlim([xMin, xMax])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("interpolation_scaling.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    interpolation_scaling()
