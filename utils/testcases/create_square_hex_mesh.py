from netCDF4 import Dataset
import numpy as np
import subprocess
import os
import argparse
from math import sqrt, cos, sin, ceil, pi
from tqdm import tqdm

#-------------------------------------------------------------------------------

def find_vertex(xV, yV, xVertex, yVertex, epsilon):

    x = np.array(xVertex)
    y = np.array(yVertex)
    d2 = (x - xV)**2 + (y - yV)**2
    matches = np.where(d2 <= epsilon**2)[0]

    return matches[0] if matches.size > 0 else None

#-------------------------------------------------------------------------------

def create_square_hex_mesh(dc,
                           lx, ly,
                           x0, y0):

    dx = (sqrt(3.0) / 2.0) * dc
    dy = dc
    a = dc / sqrt(3.0)

    nx = ceil(lx/dx)
    ny = ceil(ly/dy)

    vertexDegree = 3

    filenameOut = "grid.hex.%ix%i.nc" %(nx,ny)
    print("create: "+filenameOut)

    xCell = []
    yCell = []
    for iy in range(0, ny):
        for ix in range(0, nx):
            xCell.append(x0 + a + dx*ix)
            if (ix % 2 == 0):
                yCell.append(y0 + dc * (0.5+iy))
            else:
                yCell.append(y0 + dc * (1.0+iy))

    nCells = len(xCell)
    xCell = np.array(xCell)
    yCell = np.array(yCell)
    zCell = np.zeros(nCells)


    xVertex = []
    yVertex = []
    cellsOnVertexList = []

    epsilon = 1e-6

    iVertex = 0
    for iCell in tqdm(range(0,nCells)):
        for iVertexOnCell in range(0,6):
            angle = iVertexOnCell * ((2.0 * pi) / 6.0)
            xV = xCell[iCell] + a * cos(angle)
            yV = yCell[iCell] + a * sin(angle)

            iVertexFind = find_vertex(xV, yV, xVertex, yVertex, epsilon)
            if (iVertexFind is None):
                # vertex not found
                xVertex.append(xV)
                yVertex.append(yV)
                cellsOnVertexList.append([iCell])
                iVertex += 1
            else:
                # vertex found
                cellsOnVertexList[iVertexFind].append(iCell)

    nVertices = len(xVertex)
    xVertex = np.array(xVertex)
    yVertex = np.array(yVertex)
    zVertex = np.zeros(nVertices)

    cellsOnVertex = np.zeros((nVertices,3),dtype="i")
    for iVertex in range(0,nVertices):
        for iCellOnVertex in range(0,len(cellsOnVertexList[iVertex])):
            cellsOnVertex[iVertex,iCellOnVertex] = cellsOnVertexList[iVertex][iCellOnVertex]+1


    # create mesh file description
    filenameTmp1Out = "grid.tmp1.nc"
    fileOut = Dataset(filenameTmp1Out, "w", format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("vertexDegree", vertexDegree)

    var = fileOut.createVariable("xCell", "d", dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileOut.createVariable("yCell", "d", dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileOut.createVariable("zCell", "d", dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileOut.createVariable("xVertex", "d", dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileOut.createVariable("yVertex", "d", dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileOut.createVariable("zVertex", "d", dimensions=["nVertices"])
    var[:] = zVertex[:]

    var = fileOut.createVariable("cellsOnVertex", "d", dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileOut.on_a_sphere = "NO"

    fileOut.dc = dc
    fileOut.nx = nx
    fileOut.ny = ny
    fileOut.lx = lx
    fileOut.ly = ly
    fileOut.x0 = x0
    fileOut.y0 = y0

    fileOut.close()

    # create full mesh with the mesh converter
    MPAS_TOOLS_DIR = os.environ.get('MPAS_TOOLS_DIR')
    if (MPAS_TOOLS_DIR is None):
        raise Exception("MPAS_TOOLS_DIR environment variable nust be set")
    MpasMeshConverter = MPAS_TOOLS_DIR+"/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x"
    if (not os.path.isfile(MpasMeshConverter)):
        raise Exception("MpasMeshConverter executable must be built")

    filenameTmp2Out = filenameOut
    subprocess.run([MpasMeshConverter,filenameTmp1Out,filenameTmp2Out])

    # add attributes
    filein = Dataset(filenameOut,"a")

    filein.dc = dc
    filein.nx = nx
    filein.ny = ny
    filein.lx = lx
    filein.ly = ly
    filein.x0 = x0
    filein.y0 = y0

    filein.close()

    return filenameOut

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--dc', dest='dc', type=float, required=True)

    parser.add_argument('--lx', dest='lx', type=float, required=True)
    parser.add_argument('--ly', dest='ly', type=float, required=True)

    parser.add_argument('--x0', dest='x0', type=float, default=0.0)
    parser.add_argument('--y0', dest='y0', type=float, default=0.0)

    args = parser.parse_args()

    create_square_hex_mesh(args.dc,
                           args.lx, args.ly,
                           args.x0, args.y0)
