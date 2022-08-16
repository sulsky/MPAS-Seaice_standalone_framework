import argparse
import numpy as np
import os
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def randomize_delta_on_sphere(x, y, z, ds):

    r = np.array([x,y,z])
    rmag = np.linalg.norm(r)
    a = np.array([np.random.uniform(-1.0,1.0),
                  np.random.uniform(-1.0,1.0),
                  np.random.uniform(-1.0,1.0)])
    e1 = np.cross(r,a)
    e2 = np.cross(e1,r)
    norm1 = np.linalg.norm(e1)
    norm2 = np.linalg.norm(e2)
    e1[:] /= norm1
    e2[:] /= norm2

    ds1 = ds * np.random.uniform(-1.0,1.0)
    ds2 = ds * np.random.uniform(-1.0,1.0)

    dx = ds1 * e1[0] + ds2 * e2[0]
    dy = ds1 * e1[1] + ds2 * e2[1]
    dz = ds1 * e1[2] + ds2 * e2[2]

    return dx, dy, dz

#-------------------------------------------------------------------------------

def randomize_mesh(meshType, meshScale):

    if (meshType == "regular"):
        testName = "regular"
    elif (meshType == "random"):
        testName = "random_%f" %(meshScale)

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    meshes = [2562,10242,40962,163842]

    for mesh in meshes:

        gridnameIn = "grid.%i.nc" %(mesh)
        gridnameIntermediate = "grid.inter.%i.nc" %(mesh)
        os.system("cp %s %s" %(gridnameIn,gridnameIntermediate))

        filein = Dataset(gridnameIntermediate,"a")

        nCells = len(filein.dimensions["nCells"])
        nVertices = len(filein.dimensions["nVertices"])

        dcEdge = filein.variables["dcEdge"][:]
        avgDcEdge = np.mean(dcEdge)

        varxCell = filein.variables["xCell"]
        varyCell = filein.variables["yCell"]
        varzCell = filein.variables["zCell"]

        xCell = varxCell[:]
        yCell = varyCell[:]
        zCell = varzCell[:]

        varxVertex = filein.variables["xVertex"]
        varyVertex = filein.variables["yVertex"]
        varzVertex = filein.variables["zVertex"]

        xVertex = varxVertex[:]
        yVertex = varyVertex[:]
        zVertex = varzVertex[:]

        if (meshType == "random"):

            for iCell in range(0,nCells):
                normOld = np.linalg.norm([xCell[iCell], yCell[iCell], zCell[iCell]])
                dx, dy, dz = randomize_delta_on_sphere(xCell[iCell], yCell[iCell], zCell[iCell], avgDcEdge*meshScale)
                xCell[iCell] += dx
                yCell[iCell] += dy
                zCell[iCell] += dz
                normNew = np.linalg.norm([xCell[iCell], yCell[iCell], zCell[iCell]])
                xCell[iCell] *= (normOld/normNew)
                yCell[iCell] *= (normOld/normNew)
                zCell[iCell] *= (normOld/normNew)

            for iVertex in range(0,nVertices):
                normOld = np.linalg.norm([xVertex[iVertex], yVertex[iVertex], zVertex[iVertex]])
                dx, dy, dz = randomize_delta_on_sphere(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], avgDcEdge*meshScale)
                xVertex[iVertex] += dx
                yVertex[iVertex] += dy
                zVertex[iVertex] += dy
                normNew = np.linalg.norm([xVertex[iVertex], yVertex[iVertex], zVertex[iVertex]])
                xVertex[iVertex] *= (normOld/normNew)
                yVertex[iVertex] *= (normOld/normNew)
                zVertex[iVertex] *= (normOld/normNew)

        varxCell[:] = xCell[:]
        varyCell[:] = yCell[:]
        varzCell[:] = zCell[:]

        varxVertex[:] = xVertex[:]
        varyVertex[:] = yVertex[:]
        varzVertex[:] = zVertex[:]

        filein.close()

        gridnameOut = "grid.%s.%i.nc" %(testName,mesh)

        os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x %s %s" %(mpas_tools_dir,gridnameIntermediate,gridnameOut))

    return testName

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='meshType', default="regular", choices=['regular','random'], help='')
    parser.add_argument('-s', dest='meshScale', type=float, help='')

    args = parser.parse_args()

    randomize_mesh()
