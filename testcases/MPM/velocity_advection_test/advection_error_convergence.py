import sys

sys.path.append("../../../utils/testcases")
from log_messages import log_message

from parse_time_string import parse_time_string

from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob
import argparse

#--------------------------------------------------------

def L2_norm(numerical, analytical, nPoints, areaCell):

    norm  = 0.0
    denom = 0.0

    for iPoint in range(0,nPoints):

        norm  += areaCell[iPoint] * math.pow(numerical[iPoint] - analytical[iPoint],2)

        denom += areaCell[iPoint] * math.pow(analytical[iPoint],2)

    norm = math.sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_area(filenameMesh, filenameParticlesTemplate):

    fileMPAS = Dataset(filenameMesh, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    areaCell = fileMPAS.variables["areaCell"][:]

    fileMPAS.close()

    filenamesParticle = sorted(glob.glob(filenameParticlesTemplate))

    iceAreaCellInitial = average_particle_ice_area_to_cell(filenamesParticle[ 0], nCells)
    iceAreaCellFinal   = average_particle_ice_area_to_cell(filenamesParticle[-1], nCells)

    norm = L2_norm(iceAreaCellFinal, iceAreaCellInitial, nCells, areaCell)

    return norm
#-------------------------------------------------------------------------------

def latlon_from_xyz(x, y, z, r):

    # given xyz coordinates determine the latitude and longitude

    lon = math.atan2(y, x)
    lat = math.asin(z/r)

    return lat, lon

#--------------------------------------------------------

def posn_error(filenameTemplate):

    filenames = sorted(glob.glob(filenameTemplate))
    filein = Dataset(filenames[0],"r")
    nParticles = len(filein.dimensions["nParticles"])
    posnMP0 = filein.variables["posnMP"][0,:,:]
    filein.close()

    filein = Dataset(filenames[-1],"r")
    posnMP = filein.variables["posnMP"][0,:,:]
    xtime = filein.variables["xtime"][:,:]
    t = parse_time_string(xtime)

    r = math.sqrt(math.pow(posnMP0[0,0],2) + math.pow(posnMP0[0,1],2) + math.pow(posnMP0[0,2],2))

    normLat = 0
    normLon = 0
    for iParticle in range(0, nParticles):
        lat0, lon0 = latlon_from_xyz(posnMP0[iParticle,0], posnMP0[iParticle,1], posnMP0[iParticle,2], r)
        lat, lon = latlon_from_xyz(posnMP[iParticle,0], posnMP[iParticle,1], posnMP[iParticle,2], r)
        normLat = normLat + math.pow((lat - lat0), 2)
        normLon = normLon + math.pow((lon - (lon0 + t/r)), 2)

    normLat = math.sqrt(normLat) / nParticles
    normLon = math.sqrt(normLon) / nParticles
    return normLat, normLon

#--------------------------------------------------------

def L2_norm_particle(numerical, analytical, nParticles):

    norm  = 0.0

    for iParticle in range(0, nParticles):

            norm  = norm + math.pow(numerical[iParticle] - analytical[iParticle],2)

    norm = math.sqrt(norm) / nParticles

    return norm

#--------------------------------------------------------

def get_norm_particle(fileIC, filenameTemplate):

    filein = Dataset(fileIC,"r")
    nParticles = len(filein.dimensions["nParticles"])
    uvVelMP0 = filein.variables["uvVelMP"][:,:]
    filein.close()

    filenames = sorted(glob.glob(filenameTemplate))
    filein = Dataset(filenames[-1],"r")
    uvVelMP = filein.variables["uvVelMP"][0,:,:]
    filein.close()

    normU = L2_norm_particle(uvVelMP[:,0], uvVelMP0[:,0], nParticles)
    normV = L2_norm_particle(uvVelMP[:,1], uvVelMP0[:,1], nParticles)

    return normU, normV

#--------------------------------------------------------

def get_resolution(filename):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    # mean distance between cells
    dcEdge = fileMPAS.variables["dcEdge"][:]

    resolution = 0.0
    for iEdge in range(0,nEdges):
        resolution = resolution + dcEdge[iEdge]

    resolution = resolution / float(nEdges)

    fileMPAS.close()

    return resolution / 1000.0

#--------------------------------------------------------

def get_grid_size(filename):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    return nCells

#--------------------------------------------------------

def advection_error_convergence(runtype,
                                logFile=None):

    resolutions = [2562,10242,40962,163842]

    experiments = ["cos_lat"]

    legendLabels = ["lat", "lon", "uMP", "vMP"]

    markers = ['o','o','^','^']
    linestyles = ["-", "--", "-", "--"]
    dashes=[(1,0),(5,2.5),(1,0),(5,2.5)]

    xMin = 60
    xMax = 120

    yMin = 3e-2
    yMax = 1

    scalePos1 = 0.15e-3 / math.pow(xMax,1)
    scaleMinPos1 = math.pow(xMin,1) * scalePos1
    scaleMaxPos1 = math.pow(xMax,1) * scalePos1

    scalePos2 = 0.15e-3 / math.pow(xMax,2)
    scaleMinPos2 = math.pow(xMin,2) * scalePos2
    scaleMaxPos2 = math.pow(xMax,2) * scalePos2

    #plot
    cm = 1/2.54  # centimeters in inches
    plt.rcParams["font.family"] = "Times New Roman"
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    #positions, velocity

    fig, axes = plt.subplots(1, 1, figsize=(8*cm,6.5*cm))

    axes.loglog([xMin, xMax], [scaleMinPos1, scaleMaxPos1], linestyle=':', color='k', label="_nolegend_", lw=1)
    axes.loglog([xMin, xMax], [scaleMinPos2, scaleMaxPos2], linestyle=':', color='k', label="_nolegend_", lw=1)

    iPlot = 0
    for experiment in experiments:

         xPos = []
         yPos1 = []
         yPos2 = []
         yPos3 = []
         yPos4 = []

         for resolution in resolutions:

             fileIC = "particles_%i.nc" %(resolution)
             filename = "./output_%s_%i_%s/output.2000.nc" %(experiment,resolution,runtype)
             filenameParticlesTemplate = "./output_%s_%i_%s/particles_output*" %(experiment,resolution,runtype)

             normLat, normLon = posn_error(filenameParticlesTemplate)
             xPos.append(get_resolution(filename))
             yPos1.append(normLat)
             yPos2.append(normLon)
             normU, normV = get_norm_particle(fileIC, filenameParticlesTemplate)
             yPos3.append(normU)
             yPos4.append(normV)

         error = float('nan')
         print("error in lat: ", xPos, yPos1)
         print("error in lon: ", xPos, yPos2)
         print("error in uVel: ", xPos, yPos3)
         print("error in vVel: ", xPos, yPos4)

         axes.loglog(xPos, yPos1, marker=markers[iPlot], dashes=dashes[iPlot], color="black", markersize=5.0)
         iPlot = iPlot + 1
         axes.loglog(xPos, yPos2, marker=markers[iPlot], dashes=dashes[iPlot], color="black", markersize=5.0)
         iPlot = iPlot + 1
         axes.loglog(xPos, yPos3, marker=markers[iPlot], dashes=dashes[iPlot], color="black", markersize=5.0)
         iPlot = iPlot + 1
         axes.loglog(xPos, yPos4, marker=markers[iPlot], dashes=dashes[iPlot], color="black", markersize=5.0)


    axes.legend(legendLabels, frameon=False, loc=4, fontsize=8, handlelength=4)

    plt.minorticks_off()

    axes.set_xlabel("Grid resolution (km)")
    axes.set_ylabel(r"L2 error norm")
    axes.set_xticks([60,120,240,480])
    axes.set_xticklabels(["60","120","240","480"])
    axes.tick_params(
        axis='x',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')

    plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
    plt.savefig("advection_error_pos_convergence_%s.png" %(runtype),dpi=300)
    plt.savefig("advection_error_pos_convergence_%s.eps" %(runtype))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Plot the advection error as a function of mesh size')

    parser.add_argument('-t', required=True, dest='runtype', help='plot data for polympo or nonpolympo run')

    args = parser.parse_args()

    advection_error_convergence(args.runtype)
