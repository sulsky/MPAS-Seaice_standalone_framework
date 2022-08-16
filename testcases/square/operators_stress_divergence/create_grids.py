from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon
import argparse

#-------------------------------------------------------------------------------

def create_grid_hex(testName, meshType, meshScale, nx, ny, dc):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    Lx = float(nx) * dc
    Ly = float(ny) * dc

    gridname = "grid_%s_hex_%4.4ix%4.4i.nc" %(testName,nx,ny)

    fileout = open("namelist.input","w")

    line = "&periodic_grid\n"
    fileout.write(line)

    line = "   nx = %i,\n" %(nx+2)
    fileout.write(line)

    line = "   ny = %i,\n" %(ny+2)
    fileout.write(line)

    line = "   dc = %f,\n" %(dc)
    fileout.write(line)

    line = "   nVertLevels = 1,\n"
    fileout.write(line)

    line = "   nTracers = 1,\n"
    fileout.write(line)

    line = "   nproc = 1,\n"
    fileout.write(line)

    line = "/\n"
    fileout.write(line)

    fileout.close()

    os.system("%s/mesh_tools/periodic_hex/periodic_grid" %(mpas_tools_dir))
    os.system("python %s/mesh_tools/periodic_hex/mark_periodic_boundaries_for_culling.py -f grid.nc" %(mpas_tools_dir))
    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid.nc grid_culled.nc" %(mpas_tools_dir))

    filein = Dataset("grid_culled.nc","a")

    filein.Lx = Lx
    filein.Ly = Ly

    filein.nx = nx
    filein.ny = ny

    filein.dc = dc

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    varxCell = filein.variables["xCell"]
    varyCell = filein.variables["yCell"]

    xCell = varxCell[:]
    yCell = varyCell[:]

    varxVertex = filein.variables["xVertex"]
    varyVertex = filein.variables["yVertex"]

    xVertex = varxVertex[:]
    yVertex = varyVertex[:]

    if (meshType == "random"):

        for iCell in range(0,nCells):
            dx = np.random.uniform(-dc*meshScale,dc*meshScale)
            dy = np.random.uniform(-dc*meshScale,dc*meshScale)
            xCell[iCell] += dx
            yCell[iCell] += dy

        for iVertex in range(0,nVertices):
            dx = np.random.uniform(-dc*meshScale,dc*meshScale)
            dy = np.random.uniform(-dc*meshScale,dc*meshScale)
            xVertex[iVertex] += dx
            yVertex[iVertex] += dy

    varxCell[:] = xCell[:]
    varyCell[:] = yCell[:]

    varxVertex[:] = xVertex[:]
    varyVertex[:] = yVertex[:]

    filein.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_culled.nc %s" %(mpas_tools_dir,gridname))

#-------------------------------------------------------------------------------

def create_grid_quad(testName, meshType, meshScale, nx, ny, dc):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    Lx = float(nx) * dc
    Ly = float(ny) * dc

    gridname = "grid_%s_quad_%4.4ix%4.4i.nc" %(testName,nx,ny)

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

    if (meshType == "random"):

        for iCell in range(0,nCells):
            dx = np.random.uniform(-dc*meshScale,dc*meshScale)
            dy = np.random.uniform(-dc*meshScale,dc*meshScale)
            xCell[iCell] += dx
            yCell[iCell] += dy

        for iVertex in range(0,nVertices):
            dx = np.random.uniform(-dc*meshScale,dc*meshScale)
            dy = np.random.uniform(-dc*meshScale,dc*meshScale)
            xVertex[iVertex] += dx
            yVertex[iVertex] += dy

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

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.nc %s" %(mpas_tools_dir,gridname))

    filein = Dataset(gridname,"a")

    filein.Lx = Lx
    filein.Ly = Ly

    filein.nx = nx
    filein.ny = ny

    filein.dc = dc

    filein.close()

#-------------------------------------------------------------------------------

def create_grids(meshType, meshScale):

    if (meshType == "regular"):
        testName = "regular"
    elif (meshType == "random"):
        testName = "random_%f" %(meshScale)

    nGrid = 4
    #nGrid = 1

    # hex
    dc = 0.0125
    nx = 82
    ny = 94

    nxs = []
    nys = []
    dcs = []
    for i in range(0,nGrid):
        nxs.append(nx)
        nys.append(ny)
        dcs.append(dc)
        nx = nx*2
        ny = ny*2
        dc = dc/2

    for nx, ny, dc in zip(nxs, nys, dcs):
        print("&"*80)
        print("& Create Hex mesh: nx: ", nx, ", ny: ", ny, ", dc: ", dc)
        print("&"*80)
        create_grid_hex(testName, meshType, meshScale, nx, ny, dc)

    # quad
    dc = 0.0125
    nx = 80
    ny = 80

    nxs = []
    nys = []
    dcs = []
    for i in range(0,nGrid):
        nxs.append(nx)
        nys.append(ny)
        dcs.append(dc)
        nx = nx*2
        ny = ny*2
        dc = dc/2

    for nx, ny, dc in zip(nxs, nys, dcs):
        print("&"*80)
        print("& Create Quad mesh: nx: ", nx, ", ny: ", ny, ", dc: ", dc)
        print("&"*80)
        create_grid_quad(testName, meshType, meshScale, nx, ny, dc)

    return testName

#-------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='meshType', default="regular", choices=['regular','random'], help='')
    parser.add_argument('-s', dest='meshScale', type=float, help='')

    args = parser.parse_args()

    create_grids(args.meshType, args.meshScale)
