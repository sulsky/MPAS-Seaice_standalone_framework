from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def create_grid(gridfileOut, icfileOut):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    dc = 10000.0
    nx = 100
    ny = 100
    lx = dc * nx
    ly = dc * ny
    vertexDegree = 4

    fileGrid = Dataset("grid_in.nc","w",format="NETCDF3_CLASSIC")

    fileGrid.on_a_sphere = "NO"

    nCells = nx * ny
    nVertices = (nx+1) * (ny+1)

    fileGrid.createDimension("nCells", nCells)
    fileGrid.createDimension("nVertices", nVertices)
    fileGrid.createDimension("vertexDegree", vertexDegree)

    xCell = np.zeros(nCells)
    yCell = np.zeros(nCells)
    zCell = np.zeros(nCells)

    for j in range(0,ny):
        for i in range(0,nx):
            iCell = i + j * nx
            xCell[iCell] = (float(i) + 0.5) * dc
            yCell[iCell] = (float(j) + 0.5) * dc

    xVertex = np.zeros(nVertices)
    yVertex = np.zeros(nVertices)
    zVertex = np.zeros(nVertices)
    cellsOnVertex = np.zeros((nVertices,vertexDegree),dtype="i")

    for j in range(0,ny+1):
        for i in range(0,nx+1):
            iVertex = i + j * (nx+1)
            xVertex[iVertex] = float(i) * dc
            yVertex[iVertex] = float(j) * dc

            ic = i - 1 ; jc = j - 1
            iCell = ic + jc * nx
            if (ic >= 0 and jc >= 0):
                cellsOnVertex[iVertex,0] = iCell
            else:
                cellsOnVertex[iVertex,0] = -1

            ic = i ; jc = j - 1
            iCell = ic + jc * nx
            if (ic < nx and jc >= 0):
                cellsOnVertex[iVertex,1] = iCell
            else:
                cellsOnVertex[iVertex,1] = -1

            ic = i ; jc = j
            iCell = ic + jc * nx
            if (ic < nx and jc < ny):
                cellsOnVertex[iVertex,2] = iCell
            else:
                cellsOnVertex[iVertex,2] = -1

            ic = i - 1 ; jc = j
            iCell = ic + jc * nx
            if (ic >= 0 and jc < ny):
                cellsOnVertex[iVertex,3] = iCell
            else:
                cellsOnVertex[iVertex,3] = -1

            #print(iVertex+1,cellsOnVertex[iVertex,:]+1)

    var = fileGrid.createVariable("xCell","d",dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileGrid.createVariable("yCell","d",dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileGrid.createVariable("zCell","d",dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileGrid.createVariable("xVertex","d",dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileGrid.createVariable("yVertex","d",dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileGrid.createVariable("zVertex","d",dimensions=["nVertices"])
    var[:] = zVertex[:]

    cellsOnVertex[:] = cellsOnVertex[:] + 1
    var = fileGrid.createVariable("cellsOnVertex","i",dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileGrid.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.nc %s" %(mpas_tools_dir, gridfileOut))


    # grid
    filein = Dataset(gridfileOut,"a")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    print("%s: nCells: %i" %(gridfileOut,nCells))

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    filein.close()

    # initial v velocity
    vVelocity = np.zeros(nVertices)
    for iVertex in range(0,nVertices):
        vVelocity[iVertex] = math.sin(2.0 * math.pi * (xVertex[iVertex]/lx))
    

    # initial conditions
    fileic = Dataset(icfileOut,"w",format="NETCDF3_CLASSIC")

    fileic.createDimension("nVertices",nVertices)

    var = fileic.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = 0.0

    var = fileic.createVariable("vVelocity","d",dimensions=["nVertices"])
    var[:] = vVelocity[:]

    fileic.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_grid("grid.nc","ic.nc")
