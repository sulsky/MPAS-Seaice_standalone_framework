from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
sys.path.append("../strain")
from strain_map import get_mpas_patch_collection
from strain_map import plot_subfigure

#---------------------------------------------------------------

def strain_stress_divergence_map():

    # grid
    fileGrid = Dataset("grid.40962.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    # ic
    fileIC = Dataset("ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_wachspress_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # Weak
    fileWeak = Dataset("./output_weak_40962/output.2000.nc","r")

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("Weak:  ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()

    # mpmvar
    fileMPM = Dataset("./output_mpmvar_40962/output.2000.nc","r")

    stressDivergenceUmpmvar = fileMPM.variables["stressDivergenceU"][0,:]
    stressDivergenceVmpmvar = fileMPM.variables["stressDivergenceV"][0,:]

    stressDivergenceUmpmvarDiff = (stressDivergenceUmpmvar - stressDivergenceUAnalytical)
    stressDivergenceVmpmvarDiff = (stressDivergenceVmpmvar - stressDivergenceVAnalytical)

    print("MPMvar:  ",
          np.amin(stressDivergenceUmpmvarDiff), np.amax(stressDivergenceUmpmvarDiff),
          np.amin(stressDivergenceVmpmvarDiff), np.amax(stressDivergenceVmpmvarDiff))

    fileMPM.close()

    # mpmweak
    fileMPM = Dataset("./output_mpmweak_40962/output.2000.nc","r")

    stressDivergenceUmpmweak = fileMPM.variables["stressDivergenceU"][0,:]
    stressDivergenceVmpmweak = fileMPM.variables["stressDivergenceV"][0,:]

    stressDivergenceUmpmweakDiff = (stressDivergenceUmpmweak - stressDivergenceUAnalytical)
    stressDivergenceVmpmweakDiff = (stressDivergenceVmpmweak - stressDivergenceVAnalytical)

    print("MPMweak:  ",
          np.amin(stressDivergenceUmpmweakDiff), np.amax(stressDivergenceUmpmweakDiff),
          np.amin(stressDivergenceVmpmweakDiff), np.amax(stressDivergenceVmpmweakDiff))

    fileMPM.close()

    # mpm
    fileMPM = Dataset("./output_mpm_40962/output.2000.nc","r")

    stressDivergenceUmpm = fileMPM.variables["stressDivergenceU"][0,:]
    stressDivergenceVmpm = fileMPM.variables["stressDivergenceV"][0,:]

    stressDivergenceUmpmDiff = (stressDivergenceUmpm - stressDivergenceUAnalytical)
    stressDivergenceVmpmDiff = (stressDivergenceVmpm - stressDivergenceVAnalytical)

    print("MPM:  ",
          np.amin(stressDivergenceUmpmDiff), np.amax(stressDivergenceUmpmDiff),
          np.amin(stressDivergenceVmpmDiff), np.amax(stressDivergenceVmpmDiff))

    fileMPM.close()


    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    minVelocity = -1.0
    maxVelocity =  1.0

    minStressDiv = -20.0
    maxStressDiv =  20.0

    minStressDivDiff = -1.0
    maxStressDivDiff =  1.0

    fig, axes = plt.subplots(6, 4)

    fig.set_size_inches(7, 6.75)

    plot_subfigure(axes[0,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, uVelocity, minVelocity, maxVelocity, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$u^\prime$', '(a)', False)
    plot_subfigure(axes[0,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{u^\prime}$ Analytical', r'(b)$\times20$', False)
    plot_subfigure(axes[0,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, vVelocity, minVelocity, maxVelocity, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$v^\prime$', '(c)', False)
    plot_subfigure(axes[0,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVAnalytical, minStressDiv, maxStressDiv, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'$(\nabla \cdot \sigma)_{v^\prime}$ Analytical', r'(d)$\times20$', True)

    plot_subfigure(axes[1,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPASvar ($u^\prime$ direction)', '(e)', False)
    plot_subfigure(axes[2,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPASweak ($u^\prime$ direction)', '(e)', False)
    plot_subfigure(axes[3,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUmpmvarDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPMvar ($u^\prime$ direction)', '(q)', False)
    plot_subfigure(axes[4,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUmpmweakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPMweak ($u^\prime$ direction)', '(q)', False)
    plot_subfigure(axes[5,0], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUmpmDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM ($u^\prime$ direction)', '(q)', False)

    plot_subfigure(axes[1,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWachDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPASvar ($u^\prime$ direction)', '(f)', False)
    plot_subfigure(axes[2,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUWeakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPASweak ($u^\prime$ direction)', '(f)', False)
    plot_subfigure(axes[3,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUmpmvarDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPMvar ($u^\prime$ direction)', '(j)', False)
    plot_subfigure(axes[4,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUmpmweakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPMweak($u^\prime$ direction)', '(n)', False)
    plot_subfigure(axes[5,1], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceUmpmDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPM ($u^\prime$ direction)', '(r)', False)

    plot_subfigure(axes[1,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWachDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPASvar ($v^\prime$ direction)', '(g)', False)
    plot_subfigure(axes[2,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWeakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPASweak ($v^\prime$ direction)', '(g)', False)
    plot_subfigure(axes[3,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVmpmvarDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPMvar ($v^\prime$ direction)', '(k)', False)
    plot_subfigure(axes[4,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVmpmweakDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPMweak ($v^\prime$ direction)', '(o)', False)
    plot_subfigure(axes[5,2], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVmpmDiff, minStressDivDiff, maxStressDivDiff, -1.0, 1.0, -1.0, 1.0, \
                   False, False, r'MPM ($v^\prime$ direction)', '(s)', False)

    plot_subfigure(axes[1,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWachDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPASvar ($v^\prime$ direction)', '(h)', True)
    plot_subfigure(axes[2,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVWeakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPASweak ($v^\prime$ direction)', '(h)', True)
    plot_subfigure(axes[3,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVmpmvarDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPMvar ($v^\prime$ direction)', '(l)', True)
    plot_subfigure(axes[4,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVmpmweakDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPMweak ($v^\prime$ direction)', '(p)', True)
    plot_subfigure(axes[5,3], fig, nVertices, vertexDegreeArr, cellsOnVertex, xCell, yCell, zCell, latVertex, stressDivergenceVmpmvarDiff, minStressDivDiff, maxStressDivDiff, -0.2, 0.2, -0.2, 0.2, \
                   False, False, r'MPM ($v^\prime$ direction)', '(t)', True)

    plt.savefig("strain_stress_divergence_map.png",dpi=400)

    plt.clf()
    plt.cla()
    plt.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_stress_divergence_map()
