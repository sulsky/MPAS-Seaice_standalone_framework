from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import math, sys
import numpy as np

from matplotlib.patches import Circle
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#------------------------------------------------------------

def generate_ellipse(nEllipse):

    xEllipse = []
    yEllipse = []

    theta = math.pi / 4.0

    a = math.sqrt(2.0) * 0.5
    b = a / 2.0

    for i in range(0,nEllipse+1):

        t = (i / float(nEllipse)) * 2.0 * math.pi

        x = a * math.cos(t)
        y = b * math.sin(t)

        xp = x * math.cos(theta) - y * math.sin(theta)
        yp = x * math.sin(theta) + y * math.cos(theta)

        xp = xp - 0.5
        yp = yp - 0.5

        xEllipse.append(xp)
        yEllipse.append(yp)

    return xEllipse, yEllipse

#------------------------------------------------------------

def plot_yield_curve(axes, sig1, sig2, title, subfigureLabel, addAxisLabels=True):

    axes.scatter(sig1, sig2, s=0.08, color="black", rasterized=True)

    axes.plot((0.0,0.0), (-1.2,0.2), '--', color="grey", zorder=2, linewidth=0.5)
    axes.plot((-1.2,0.2), (0.0,0.0), '--', color="grey", zorder=3, linewidth=0.5)

    xEllipse, yEllipse = generate_ellipse(100)

    axes.tick_params('both', length=3)

    axes.plot(xEllipse,yEllipse, ':', color="grey", zorder=4, linewidth=0.5)

    axes.set_xlim([-1.2,0.2])
    axes.set_ylim([-1.2,0.2])

    if (not addAxisLabels):
        axes.get_xaxis().set_ticklabels([])
        axes.get_yaxis().set_ticklabels([])

    axes.set_aspect('equal')

    axes.set_title(subfigureLabel+" "+title,loc='left')

    if (addAxisLabels):
        axes.set_xlabel(r"$\sigma_1$")
        axes.set_ylabel(r"$\sigma_2$")

#------------------------------------------------------------

def get_cice_stresses(filename):

    fileCICE = Dataset(filename, "r")

    sig1in = fileCICE.variables["sig1"][0,:,:]
    sig2in = fileCICE.variables["sig2"][0,:,:]

    fileCICE.close()

    sig1 = []
    sig2 = []

    for j in range(0,84):
        for i in range(0,84):
            if (sig1in[j,i] < 1e20):
                sig1.append(sig1in[j,i])
                sig2.append(sig2in[j,i])

    return sig1, sig2

#------------------------------------------------------------

def get_mpas_variational_stresses(filename):

    iTime = 3

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    nEdgesOnCell = fileMPAS.variables["nEdgesOnCell"][:]

    principalStress1 = fileMPAS.variables["principalStress1Var"][iTime,0:nCells-1,:]
    principalStress2 = fileMPAS.variables["principalStress2Var"][iTime,0:nCells-1,:]

    sig1 = []
    sig2 = []

    for iCell in range(0,nCells-1):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            if (principalStress1[iCell,iVertexOnCell] < 1e20):
                sig1.append(principalStress1[iCell,iVertexOnCell])
                sig2.append(principalStress2[iCell,iVertexOnCell])

    fileMPAS.close()

    return sig1, sig2

#------------------------------------------------------------

def get_mpas_weak_stresses(filename):

    iTime = 3

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    principalStress1 = fileMPAS.variables["principalStress1Weak"][iTime,0:nCells-1]
    principalStress2 = fileMPAS.variables["principalStress2Weak"][iTime,0:nCells-1]

    sig1 = []
    sig2 = []

    for iCell in range(0,nCells-1):
        if (principalStress1[iCell] < 1e20):
            sig1.append(principalStress1[iCell])
            sig2.append(principalStress2[iCell])

    fileMPAS.close()

    return sig1, sig2

#------------------------------------------------------------

def get_mpas_mpm_stresses(filename):

    fileMPAS = Dataset(filename, "r")

    nParticles = len(fileMPAS.dimensions["nParticles"])

    principalStress1 = fileMPAS.variables["principalStress1MP"][:,:]
    principalStress2 = fileMPAS.variables["principalStress2MP"][:,:]

    sig1 = []
    sig2 = []

    for iParticle in range(0,nParticles):
        if (principalStress1[0,iParticle] < 1e20):
            sig1.append(principalStress1[0,iParticle])
            sig2.append(principalStress2[0,iParticle])

    fileMPAS.close()

    return sig1, sig2

#------------------------------------------------------------

cm = 1/2.54  # centimeters in inches
plt.rc('font', family="Times New Roman")
plt.rc('mathtext',fontset="stix")
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


sig1MPASQuadWach, sig2MPASQuadWach = get_mpas_variational_stresses("./output_quad_wachspress_0080x0080_120/output.2000.nc")
sig1MPASQuadWeak, sig2MPASQuadWeak = get_mpas_weak_stresses("./output_quad_weak_0080x0080_120/output.2000.nc")

sig1MPASQuadMPM, sig2MPASQuadMPM = get_mpas_mpm_stresses("./output_quad_mpmmpm_0080x0080_120/particles_output.2000-01-01_04.00.00.nc")

sig1MPASQuadMPMWeak, sig2MPASQuadMPMWeak = get_mpas_mpm_stresses("./output_quad_mpmweak_0080x0080_120/particles_output.2000-01-01_04.00.00.nc")

sig1MPASQuadMPMVar, sig2MPASQuadMPMVar= get_mpas_mpm_stresses("./output_quad_mpmvar_0080x0080_120/particles_output.2000-01-01_04.00.00.nc")

sig1MPASHexWach,  sig2MPASHexWach = get_mpas_variational_stresses("./output_hex_wachspress_0082x0094_120/output.2000.nc")
sig1MPASHexWeak,  sig2MPASHexWeak = get_mpas_weak_stresses("./output_hex_weak_0082x0094_120/output.2000.nc")

sig1MPASHexMPM, sig2MPASHexMPM = get_mpas_mpm_stresses("./output_hex_mpmmpm_0082x0094_120/particles_output.2000-01-01_04.00.00.nc")

sig1MPASHexMPMWeak, sig2MPASHexMPMWeak = get_mpas_mpm_stresses("./output_hex_mpmweak_0082x0094_120/particles_output.2000-01-01_04.00.00.nc")

sig1MPASHexMPMVar, sig2MPASHexMPMVar = get_mpas_mpm_stresses("./output_hex_mpmvar_0082x0094_120/particles_output.2000-01-01_04.00.00.nc")

fig, axes = plt.subplots(2, 5, figsize=(15*cm,7.4*cm))

plot_yield_curve(axes[0,0], sig1MPASQuadWach, sig2MPASQuadWach, "Quad Wachs. Var", "(a)", False)
plot_yield_curve(axes[0,1], sig1MPASQuadMPMVar,  sig2MPASQuadMPMVar,  "Quad MPM Var",    "(b)", False)
plot_yield_curve(axes[0,2], sig1MPASQuadWeak, sig2MPASQuadWeak, "Quad Wachs. Weak",   "(c)", False)
plot_yield_curve(axes[0,3], sig1MPASQuadMPMWeak, sig2MPASQuadMPMWeak, "Quad MPM Weak",   "(d)", False)
plot_yield_curve(axes[0,4], sig1MPASQuadMPM,  sig2MPASQuadMPM,  "Quad MPM MPM",    "(e)", False)

plot_yield_curve(axes[1,0], sig1MPASHexWach,  sig2MPASHexWach,  "Hex Wachs. Var",  "(f)", False)
plot_yield_curve(axes[1,1], sig1MPASHexMPMVar,   sig2MPASHexMPMVar,   "Hex MPM Var",     "(g)", False)
plot_yield_curve(axes[1,2], sig1MPASHexWeak,  sig2MPASHexWeak,  "Hex Wachs. Weak",    "(h)", False)
plot_yield_curve(axes[1,3], sig1MPASHexMPMWeak,  sig2MPASHexMPMWeak,  "Hex MPM Weak",    "(i)", False)
plot_yield_curve(axes[1,4], sig1MPASHexMPM,   sig2MPASHexMPM,   "Hex MPM MPM",     "(j)", False)

plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
plt.savefig("principal_stresses.eps",dpi=300)
plt.savefig("principal_stresses.png",dpi=300)
plt.savefig("principal_stresses.pdf",dpi=300)
