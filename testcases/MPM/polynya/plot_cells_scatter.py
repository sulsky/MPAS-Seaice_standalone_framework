from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def plot_cells_scatter():

    arrayName = "iceVolumeCell"

    ymin = 200000.0
    ymax = 300000.0

    filein = Dataset("grid.nc","r")

    nCells = len(filein.dimensions["nCells"])

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    filein.close()

    indices = np.where((yCell >= ymin) & (yCell <= ymax))[0]

    filein = Dataset("./output_cells/output.2000.nc","r")

    nTimes = len(filein.dimensions["Time"])

    arrayCell = filein.variables[arrayName][:]

    filein.close()

    fig, axis = plt.subplots()

    cmap = plt.cm.jet
    colors = cmap(np.linspace(0, 1, nTimes))

    for iTime in range(0,nTimes):
        axis.scatter(xCell[indices], arrayCell[iTime,indices], color=colors[iTime], s=1, rasterized=True)

    axis.set_xlim((0,1000000))
    axis.set_ylim((0,0.6))
    plt.tight_layout()
    plt.savefig("cells_scatter.png",dpi=300)
    plt.savefig("cells_scatter.pdf",dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_cells_scatter()
