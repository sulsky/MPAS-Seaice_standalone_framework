import sys
from math import sqrt, pi, sin, cos, asin, acos
import numpy as np
import argparse
from netCDF4 import Dataset
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def mpas_arc_length(ax, ay, az, bx, by, bz):

    cx = bx - ax
    cy = by - ay
    cz = bz - az

    #      r = ax*ax + ay*ay + az*az
    #      c = cx*cx + cy*cy + cz*cz
    #
    #      arc_length = sqrt(r) * acos(1.0 - c/(2.0*r))

    r = sqrt(ax*ax + ay*ay + az*az)
    c = sqrt(cx*cx + cy*cy + cz*cz)
    #      arc_length = sqrt(r) * 2.0 * asin(c/(2.0*r))
    return r * 2.0 * asin(c/(2.0*r))

#-------------------------------------------------------------------------------

def mpas_sphere_angle(ax, ay, az, bx, by, bz, cx, cy, cz):

    a = acos(max(min(bx*cx + by*cy + bz*cz,1.0),-1.0))      # Eqn. (3)
    b = acos(max(min(ax*cx + ay*cy + az*cz,1.0),-1.0))      # Eqn. (2)
    c = acos(max(min(ax*bx + ay*by + az*bz,1.0),-1.0))      # Eqn. (1)

    ABx = bx - ax
    ABy = by - ay
    ABz = bz - az

    ACx = cx - ax
    ACy = cy - ay
    ACz = cz - az

    Dx =   (ABy * ACz) - (ABz * ACy)
    Dy = -((ABx * ACz) - (ABz * ACx))
    Dz =   (ABx * ACy) - (ABy * ACx)

    s = 0.5*(a + b + c)
    #      sin_angle = sqrt((sin(s-b)*sin(s-c))/(sin(b)*sin(c)))   # Eqn. (28)
    sin_angle = sqrt(min(1.0,max(0.0,(sin(s-b)*sin(s-c))/(sin(b)*sin(c)))))   # Eqn. (28)

    if ((Dx*ax + Dy*ay + Dz*az) >= 0.0):
        return  2.0 * asin(max(min(sin_angle,1.0),-1.0))
    else:
        return -2.0 * asin(max(min(sin_angle,1.0),-1.0))

#-------------------------------------------------------------------------------

def mpas_rotate_about_vector(x, y, z, theta, a, b, c, u, v, w):

    vw2 = v*v + w*w
    uw2 = u*u + w*w
    uv2 = u*u + v*v
    m = sqrt(u*u + v*v + w*w)

    xp = (a*vw2 + u*(-b*v-c*w+u*x+v*y+w*z) + ((x-a)*vw2+u*(b*v+c*w-v*y-w*z))*cos(theta) + m*(-c*v+b*w-w*y+v*z)*sin(theta))/(m*m)
    yp = (b*uw2 + v*(-a*u-c*w+u*x+v*y+w*z) + ((y-b)*uw2+v*(a*u+c*w-u*x-w*z))*cos(theta) + m*( c*u-a*w+w*x-u*z)*sin(theta))/(m*m)
    zp = (c*uv2 + w*(-a*u-b*v+u*x+v*y+w*z) + ((z-c)*uv2+w*(a*u+b*v-u*x-v*y))*cos(theta) + m*(-b*u+a*v-v*x+u*y)*sin(theta))/(m*m)

    return xp, yp, zp

#-------------------------------------------------------------------------------

def mpas_mirror_point(xPoint,
                      yPoint,
                      zPoint,
                      xA,
                      yA,
                      zA,
                      xB,
                      yB,
                      zB):

    #
    # Find the spherical angle between arcs AP and AB (where P is the point to be reflected)
    #
    alpha = mpas_sphere_angle(xA, yA, zA, xPoint, yPoint, zPoint, xB, yB, zB)

    #
    # Rotate the point to be reflected by twice alpha about the vector from the origin to A
    #
    xMirror, yMirror, zMirror = mpas_rotate_about_vector(xPoint,
                                                         yPoint,
                                                         zPoint,
                                                         2.0 * alpha,
                                                         0.0,
                                                         0.0,
                                                         0.0,
                                                         xA,
                                                         yA,
                                                         zA)

    return xMirror, yMirror, zMirror

#-------------------------------------------------------------------------------

def mpas_in_cell(xPoint,
                 yPoint,
                 zPoint,
                 xCell,
                 yCell,
                 zCell,
                 nEdgesOnCell,
                 verticesOnCell,
                 xVertex,
                 yVertex,
                 zVertex):

    radius = sqrt(xCell * xCell + yCell * yCell + zCell * zCell)
    radius_inv = 1.0 / radius

    inDist = mpas_arc_length(xPoint, yPoint, zPoint, xCell, yCell, zCell)

    mpas_in_cell = True

    for i in range(0,nEdgesOnCell):
        vtx1 = verticesOnCell[i]
        vtx2 = verticesOnCell[(i+1) % nEdgesOnCell]

        xNeighbor, yNeighbor, zNeighbor = mpas_mirror_point(xCell*radius_inv,
                                                            yCell*radius_inv,
                                                            zCell*radius_inv,
                                                            xVertex[vtx1]*radius_inv,
                                                            yVertex[vtx1]*radius_inv,
                                                            zVertex[vtx1]*radius_inv,
                                                            xVertex[vtx2]*radius_inv,
                                                            yVertex[vtx2]*radius_inv,
                                                            zVertex[vtx2]*radius_inv)

        xNeighbor = xNeighbor * radius
        yNeighbor = yNeighbor * radius
        zNeighbor = zNeighbor * radius

        outDist = mpas_arc_length(xPoint, yPoint, zPoint, xNeighbor, yNeighbor, zNeighbor)

        if (outDist < inDist):
            mpas_in_cell = False
            return mpas_in_cell

    return mpas_in_cell

#-------------------------------------------------------------------------------

def place_particles(posnMP,
                    latlon,
                    iCellMP,
                    creationIndexMP,
                    iCell,
                    nParticlesPerCellDesired,
                    nParticlesPerCellActual,
                    posnInitType,
                    on_a_sphere,
                    earthRadius,
                    nEdgesOnCell,
                    verticesOnCell,
                    latVertex,
                    lonVertex,
                    xVertex,
                    yVertex,
                    zVertex,
                    xCell,
                    yCell,
                    zCell):

    if (posnInitType.strip() == 'even'):

        # place points evenly in circumscribing rectangle
        # WARNING: for 'even' to work, numberToPlace must be a perfect square
        np = int(sqrt(float(nParticlesPerCellDesired)))
        #if ( (np*np) \= initNumMPPerCell(iCell)) then
        #   print*,'error in placing material points'
        #endif

        if (on_a_sphere):

            # determine rectangle
            latmin =  sys.float_info.max
            latmax = -sys.float_info.max
            lonmin =  sys.float_info.max
            lonmax = -sys.float_info.max
            for iVertexOnCell in range(0, nEdgesOnCell):
                iVertex = verticesOnCell[iVertexOnCell]
                latmin = min(latmin,latVertex[iVertex])
                latmax = max(latmax,latVertex[iVertex])
                lonmin = min(lonmin,lonVertex[iVertex])
                lonmax = max(lonmax,lonVertex[iVertex])

            dx = (latmax - latmin) / float(np) / 2.0
            dy = (lonmax - lonmin) / float(np) / 2.0
            # special case for lon near -pi and pi
            if (lonmin * lonmax < 0.0 and lonmax > pi / 2.0):
                dy = (lonmax - lonmin - 2.0 * pi) / float(np) / 2.0


            k = 0
            for i in range(0, np):
                for j in range(0, np):

                    lat = latmin + dx * float(2 * i + 1)
                    lon = lonmin + dy * float(2 * j + 1)
                    if (lon < -pi): lon = 2.0 * pi + lon

                    x = earthRadius * cos(lon) * cos(lat)
                    y = earthRadius * sin(lon) * cos(lat)
                    z = earthRadius * sin(lat)

                    if (mpas_in_cell(x,
                                     y,
                                     z,
                                     xCell,
                                     yCell,
                                     zCell,
                                     nEdgesOnCell,
                                     verticesOnCell[:],
                                     xVertex,
                                     yVertex,
                                     zVertex)):

                        k = k + 1
                        posnMP.append([x, y, z])
                        latlon.append([lat, lon])
                        iCellMP.append(iCell+1)
                        creationIndexMP.append(k)

            nParticlesPerCellActual = k

        else: # not on a sphere

            # determine rectangle
            xmin =  sys.float_info.max
            xmax = -sys.float_info.max
            ymin =  sys.float_info.max
            ymax = -sys.float_info.max
            for iVertexOnCell in range(0, nEdgesOnCell):
                iVertex = verticesOnCell[iVertexOnCell]
                xmin = min(xmin, xVertex[iVertex])
                xmax = max(xmax, xVertex[iVertex])
                ymin = min(ymin, yVertex[iVertex])
                ymax = max(ymax, yVertex[iVertex])

            dx = (xmax-xmin) / float(np) / 2.0
            dy = (ymax-ymin) / float(np) / 2.0

            k = 0
            for i in range(0, np):
                for j in range(0, np):

                    x = xmin + dx * float(2 * i + 1)
                    y = ymin + dy * float(2 * j + 1)
                    z = 0.0

                    if (mpas_in_cell(x,
                                     y,
                                     z,
                                     xCell,
                                     yCell,
                                     zCell,
                                     nEdgesOnCell,
                                     verticesOnCell[:],
                                     xVertex,
                                     yVertex,
                                     zVertex)):

                        k = k + 1
                        posnMP.append([x, y, z])
                        latlon.append([0.0, 0.0])
                        iCellMP.append(iCell+1)
                        creationIndexMP.append(k)

            nParticlesPerCellActual = k

    elif (posnInitType.strip() == 'poisson'):
        # place points using a poisson distribution
        raise Exception("Poisson position init type not implemented")
    elif (posnInitType.strip() == 'random'):
        # place points using a random distribution
        raise Exception("Random position init type not implemented")
    else:
        raise Exception("Invalid particle initType")

    return posnMP, latlon, iCellMP, creationIndexMP, nParticlesPerCellActual

#-------------------------------------------------------------------------------

def initial_particle_positions(filenameMesh,
                               filenameIceConc,
                               filenameOut,
                               particleInitType,
                               particleInitNumber,
                               particlePositionInitType,
                               earthRadius,
                               iceConcTimeIndex=-1):

    # mesh info
    fileMesh = Dataset(filenameMesh,"r")

    on_a_sphere = fileMesh.on_a_sphere

    nCells = len(fileMesh.dimensions["nCells"])

    nEdgesOnCell = fileMesh.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMesh.variables["verticesOnCell"][:]

    areaCell = fileMesh.variables["areaCell"][:]
    latVertex = fileMesh.variables["latVertex"][:]
    lonVertex = fileMesh.variables["lonVertex"][:]
    xVertex = fileMesh.variables["xVertex"][:]
    yVertex = fileMesh.variables["yVertex"][:]
    zVertex = fileMesh.variables["zVertex"][:]
    xCell = fileMesh.variables["xCell"][:]
    yCell = fileMesh.variables["yCell"][:]
    zCell = fileMesh.variables["zCell"][:]

    fileMesh.close()

    verticesOnCell[:] -= 1

    # ice concentration
    fileIceConc = Dataset(filenameIceConc,"r")

    nCellsConc = len(fileIceConc.dimensions["nCells"])

    nDims = fileIceConc.variables["iceAreaCell"].ndim
    if (nDims == 1):
        iceAreaCell = fileIceConc.variables["iceAreaCell"][:]
    elif (nDims == 2):
        iceAreaCell = fileIceConc.variables["iceAreaCell"][iceConcTimeIndex,:]
    else:
        raise Exception("Ice concentration variable has too many dimensions")

    fileIceConc.close()

    if (nCellsConc != nCells):
        raise Exception("Inconsistent cells dimension: %i %i" %(nCells,nCellsConc))

    # average cell size
    averageCellSize = np.mean(areaCell)


    # particle positions
    initNumMPPerCell = np.zeros(nCells,dtype="i")

    posnMP = []
    posnLatLonGeoMP = []
    cellIDCreationMP = []
    creationIndexMP = []

    for iCell in range(0, nCells):

        if (iceAreaCell[iCell] > 0.0):

            # TODO add some checks for minimum sizes, possibly reduce
            # number of added material points to consolidate
            if (particleInitType.strip() == 'number'):
                nParticlesPerCellDesired = particleInitNumber
            elif (particleInitType.strip() == 'area'):
                nParticlesPerCellDesired = particleInitNumber * int(areaCell[iCell] / averageCellSize)
                nParticlesPerCellDesired = max(1, nParticlesPerCellDesired)
            else:
                # error not a valid option
                raise Exception("Invalid particle_init_type")


            posnMP, latlon, iCellMP, creationIndexMP, initNumMPPerCell[iCell] = place_particles(
                posnMP,
                posnLatLonGeoMP,
                cellIDCreationMP,
                creationIndexMP,
                iCell,
                nParticlesPerCellDesired,
                initNumMPPerCell[iCell],
                particlePositionInitType,
                on_a_sphere,
                earthRadius,
                nEdgesOnCell[iCell],
                verticesOnCell[iCell,:],
                latVertex,
                lonVertex,
                xVertex,
                yVertex,
                zVertex,
                xCell[iCell],
                yCell[iCell],
                zCell[iCell])

        else:

            initNumMPPerCell[iCell] = 0

    nParticles = len(posnMP)
    posnMP = np.array(posnMP)
    posnLatLonGeoMP = np.array(posnLatLonGeoMP)
    cellIDCreationMP = np.array(cellIDCreationMP)
    creationIndexMP = np.array(creationIndexMP)

    # output
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nParticles",nParticles)
    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("THREE", 3)
    fileOut.createDimension("TWO", 2)

    var = fileOut.createVariable("posnMP","d",dimensions=["nParticles","THREE"])
    var[:] = posnMP[:]

    var = fileOut.createVariable("posnLatLonGeoMP","d",dimensions=["nParticles","TWO"])
    var[:] = posnLatLonGeoMP[:]

    var = fileOut.createVariable("cellIDCreationMP","i",dimensions=["nParticles"])
    var[:] = cellIDCreationMP[:]

    var = fileOut.createVariable("creationIndexMP","i",dimensions=["nParticles"])
    var[:] = creationIndexMP[:]

    var = fileOut.createVariable("initNumMPPerCell","i",dimensions=["nCells"])
    var[:] = initNumMPPerCell[:]

    fileOut.close()


    # plot
    fig = plt.figure(figsize=(10,10))
    axis = fig.add_subplot(projection='3d')
    axis.scatter(posnMP[:,0],posnMP[:,1],posnMP[:,2])
    axis.set_xlim((np.amin(xCell)*1.1,np.amax(xCell)*1.1))
    axis.set_ylim((np.amin(yCell)*1.1,np.amax(yCell)*1.1))
    axis.set_zlim((np.amin(zCell)*1.1,np.amax(zCell)*1.1))
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    plt.tight_layout()
    plt.savefig("particles.png")

#-------------------------------------------------------------------------------

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser()

    parser.add_argument('-m', dest="filenameMesh", required=True)
    parser.add_argument('-c', dest="filenameIceConc", required=True)
    parser.add_argument('-o', dest="filenameOut", required=True)
    parser.add_argument('-t', dest="particleInitType", required=True, choices=["number","area"])
    parser.add_argument('-n', dest="particleInitNumber", required=True, type=int)
    parser.add_argument('-p', dest="particlePositionInitType", required=True, choices=["even","poisson","random"])
    parser.add_argument('-r', dest="earthRadius", type=float, default=6371229.0)
    parser.add_argument('--concIndex', dest="iceConcTimeIndex", default=-1)

    args = parser.parse_args()

    initial_particle_positions(args.filenameMesh,
                               args.filenameIceConc,
                               args.filenameOut,
                               args.particleInitType,
                               args.particleInitNumber,
                               args.particlePositionInitType,
                               args.earthRadius,
                               args.iceConcTimeIndex)
