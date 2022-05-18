from netCDF4 import Dataset
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.linalg import lu_factor, lu_solve

#--------------------------------------------------------
# General
#--------------------------------------------------------

def local_eastern_and_northern_unit_vectors(xPoint, yPoint, zPoint):

    unitVectorEast = np.zeros(3)
    unitVectorNorth = np.zeros(3)

    # determine unit vector at cell centre pointing in U directions - local eastwards
    unitVectorEast[0] = -yPoint
    unitVectorEast[1] =  xPoint
    unitVectorEast[2] =  0.0

    vectorMagnitude = math.sqrt(math.pow(unitVectorEast[0],2) +
                                math.pow(unitVectorEast[1],2) +
                                math.pow(unitVectorEast[2],2))

    unitVectorEast[0] /= vectorMagnitude
    unitVectorEast[1] /= vectorMagnitude
    unitVectorEast[2] /= vectorMagnitude

    # determine unit vector at cell centre pointing in V directions - local northwards
    if (zPoint != 0.0):

       unitVectorNorth[0] = -xPoint
       unitVectorNorth[1] = -yPoint
       unitVectorNorth[2] = (math.pow(xPoint,2) + math.pow(yPoint,2)) / zPoint

       vectorMagnitude = math.sqrt(math.pow(unitVectorNorth[0],2) +
                                   math.pow(unitVectorNorth[1],2) +
                                   math.pow(unitVectorNorth[2],2))

       unitVectorNorth[0] /= vectorMagnitude
       unitVectorNorth[1] /= vectorMagnitude
       unitVectorNorth[2] /= vectorMagnitude

       if (zPoint < 0.0):

          # have calculated due south vector if in southern hemisphere
          unitVectorNorth[0] = -unitVectorNorth[0]
          unitVectorNorth[1] = -unitVectorNorth[1]
          unitVectorNorth[2] = -unitVectorNorth[2]

    else:

       unitVectorNorth[0] = 0.0
       unitVectorNorth[1] = 0.0
       unitVectorNorth[2] = 1.0

    return unitVectorEast, unitVectorNorth

#--------------------------------------------------------

def dot_product_3space(a, b):

    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

#--------------------------------------------------------

def project_3D_vector_onto_local_2D(normalVector3D,
                                    xPoint,
                                    yPoint,
                                    zPoint):

    normalVector2D = np.zeros(2)

    # calculate the local eastern and northern unit vectors at the point
    unitVectorEast, unitVectorNorth = local_eastern_and_northern_unit_vectors(xPoint, yPoint, zPoint)

    # project normal 3D vector onto U 2D direction - eastwards
    normalVector2D[0] = dot_product_3space(normalVector3D, unitVectorEast)

    # project normal 3D vector onto V 2D direction - northwards
    normalVector2D[1] = dot_product_3space(normalVector3D, unitVectorNorth)

    return normalVector2D

#--------------------------------------------------------

def grid_rotation_forward(x, y, z, rotateCartesianGrid):

    if (rotateCartesianGrid):

       xp = -z
       yp = y
       zp = x

    else:

       xp = x
       yp = y
       zp = z

    return xp, yp, zp

#--------------------------------------------------------

def grid_rotation_backward(vxp, vyp, vzp, rotateCartesianGrid):

    if (rotateCartesianGrid):

       vx =  vzp
       vy =  vyp
       vz = -vxp

    else:

       vx = vxp
       vy = vyp
       vz = vzp

    return vx, vy, vz

#--------------------------------------------------------

def local_coords(nCells,
                 maxEdges,
                 nEdgesOnCell,
                 verticesOnCell,
                 xCell,
                 yCell,
                 zCell,
                 xVertex,
                 yVertex,
                 zVertex,
                 rotateCartesianGrid):

    xLocal = np.zeros((maxEdges,nCells))
    yLocal = np.zeros((maxEdges,nCells))

    for iCell in range(0,nCells):

       xCellRotated, yCellRotated, zCellRotated = grid_rotation_forward(xCell[iCell], yCell[iCell], zCell[iCell], rotateCartesianGrid)

       for iVertexOnCell in range(0, nEdgesOnCell[iCell]):

          iVertex = verticesOnCell[iCell,iVertexOnCell]

          normalVector3D = grid_rotation_forward(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], rotateCartesianGrid)

          normalVector2D = project_3D_vector_onto_local_2D(normalVector3D, xCellRotated, yCellRotated, zCellRotated)

          xLocal[iVertexOnCell,iCell] = normalVector2D[0]
          yLocal[iVertexOnCell,iCell] = normalVector2D[1]

    return xLocal, yLocal

#--------------------------------------------------------

def local_to_global_coordinates(x3DC, y3DC, z3DC,
                                unitVectorEast,
                                unitVectorNorth,
                                x2D, y2D):

    x3D = x3DC + x2D * unitVectorEast[0] + y2D * unitVectorNorth[0]
    y3D = y3DC + x2D * unitVectorEast[1] + y2D * unitVectorNorth[1]
    z3D = z3DC + x2D * unitVectorEast[2] + y2D * unitVectorNorth[2]

    return x3D, y3D, z3D

#--------------------------------------------------------

def latlon_from_xyz(x, y, z, r):

    lon = math.atan2(y, x)
    lat = math.asin(z/r)

    return lat, lon

#--------------------------------------------------------

def wrapped_index(i,n):

    j = i
    if (j < 0): j += n
    if (j > n-1): j -= n

    return j

#--------------------------------------------------------
# Wachspress
#--------------------------------------------------------

def wachspress_indexes(nEdgesOnCell):

    nEdgesOnCellSubset = np.zeros(nEdgesOnCell,dtype="i")
    vertexIndexSubset = np.zeros((nEdgesOnCell,nEdgesOnCell),dtype="i")

    for jVertex in range(0,nEdgesOnCell):

        i1 = jVertex
        i2 = wrapped_index(jVertex + 1, nEdgesOnCell)

        nEdgesOnCellSubset[jVertex] = 0

        for kVertex in range(0, nEdgesOnCell):

            if (kVertex != i1 and kVertex != i2):
                nEdgesOnCellSubset[jVertex] = nEdgesOnCellSubset[jVertex] + 1
                vertexIndexSubset[jVertex,nEdgesOnCellSubset[jVertex]-1] = kVertex

    return nEdgesOnCellSubset, vertexIndexSubset

#--------------------------------------------------------

def calc_wachspress_coefficients(nCells,
                                 maxEdges,
                                 nEdgesOnCell,
                                 xLocal,
                                 yLocal):

    wachspressKappa = np.zeros((maxEdges,maxEdges,nCells))
    wachspressA = np.zeros((maxEdges,nCells))
    wachspressB = np.zeros((maxEdges,nCells))

    # loop over cells
    for iCell in range(0,nCells):

        # loop over vertices
        for iVertex in range(0,nEdgesOnCell[iCell]):

            # end points of line segment
            i1 = iVertex - 1
            i2 = iVertex
            if (i1 < 0):
                i1 = i1 + nEdgesOnCell[iCell]

            # solve for the line segment equation
            wachspressA[iVertex, iCell] = (yLocal[i2,iCell] - yLocal[i1,iCell]) / (xLocal[i1,iCell] * yLocal[i2,iCell] - xLocal[i2,iCell] * yLocal[i1,iCell])
            wachspressB[iVertex, iCell] = (xLocal[i1,iCell] - xLocal[i2,iCell]) / (xLocal[i1,iCell] * yLocal[i2,iCell] - xLocal[i2,iCell] * yLocal[i1,iCell])

        # loop over vertices
        for iVertex in range(0,nEdgesOnCell[iCell]):

            # determine kappa
            wachspressKappa[0,iVertex,iCell] = 1.0

            for jVertex in range(1, nEdgesOnCell[iCell]):

                # previous, this and next vertex
                i0 = jVertex - 1
                i1 = jVertex
                i2 = jVertex + 1
                if (i2 >= nEdgesOnCell[iCell]):
                    i2 = i2 - nEdgesOnCell[iCell]

                wachspressKappa[jVertex,iVertex,iCell] = wachspressKappa[jVertex-1,iVertex,iCell] * \
                    (wachspressA[i2,iCell] * (xLocal[i0,iCell] - xLocal[i1,iCell]) + \
                     wachspressB[i2,iCell] * (yLocal[i0,iCell] - yLocal[i1,iCell])) / \
                    (wachspressA[i0,iCell] * (xLocal[i1,iCell] - xLocal[i0,iCell]) + \
                     wachspressB[i0,iCell] * (yLocal[i1,iCell] - yLocal[i0,iCell]))

    return wachspressKappa, wachspressA, wachspressB

#--------------------------------------------------------

def wachspress_edge_equation(x, y, wachspressA, wachspressB):

    edgeEquation = 1.0 - wachspressA * x - wachspressB * y

    return edgeEquation

#--------------------------------------------------------

def wachspress_numerator(nEdgesOnCell,
                         jVertex,
                         iVertex,
                         x,
                         y,
                         wachspressKappa,
                         wachspressA,
                         wachspressB,
                         nEdgesOnCellSubset,
                         vertexIndexSubset):

    numerator = 1.0

    for kVertex in range(0,nEdgesOnCellSubset[jVertex]):

        edgeEquation = wachspress_edge_equation(x,
                                                y,
                                                wachspressA[vertexIndexSubset[jVertex,kVertex]],
                                                wachspressB[vertexIndexSubset[jVertex,kVertex]])

        numerator = numerator * edgeEquation

    numerator = numerator * wachspressKappa[jVertex,iVertex]

    return numerator

#--------------------------------------------------------

def wachspress_basis_function(nEdgesOnCell,
                              iVertex,
                              x,
                              y,
                              wachspressKappa,
                              wachspressA,
                              wachspressB,
                              nEdgesOnCellSubset,
                              vertexIndexSubset):

    numerator = np.zeros(nEdgesOnCell)
    denominator = 0.0

    for jVertex in range(0,nEdgesOnCell):

        numerator[jVertex] = wachspress_numerator(nEdgesOnCell,
                                                  jVertex,
                                                  iVertex,
                                                  x,
                                                  y,
                                                  wachspressKappa,
                                                  wachspressA,
                                                  wachspressB,
                                                  nEdgesOnCellSubset,
                                                  vertexIndexSubset)

        denominator = denominator + numerator[jVertex]

    wachpress = numerator[iVertex] / denominator

    return wachpress

#--------------------------------------------------------

def wachspress_edge_equation_array(x, y, wachspressA, wachspressB):

    edgeEquation = np.zeros(np.shape(x)[0])

    edgeEquation[:] = 1.0 - wachspressA * x[:] - wachspressB * y[:]

    return edgeEquation

#--------------------------------------------------------

def wachspress_numerator_array(nEdgesOnCell,
                               jVertex,
                               iVertex,
                               x,
                               y,
                               wachspressKappa,
                               wachspressA,
                               wachspressB,
                               nEdgesOnCellSubset,
                               vertexIndexSubset):

    numerator = np.ones(np.shape(x)[0])

    for kVertex in range(0,nEdgesOnCellSubset[jVertex]):

        edgeEquation = wachspress_edge_equation(x,
                                                y,
                                                wachspressA[vertexIndexSubset[jVertex,kVertex]],
                                                wachspressB[vertexIndexSubset[jVertex,kVertex]])

        numerator[:] = numerator[:] * edgeEquation[:]

    numerator[:] = numerator[:] * wachspressKappa[jVertex,iVertex]

    return numerator

#--------------------------------------------------------

def wachspress_basis_function_array(nEdgesOnCell,
                                    iVertex,
                                    x,
                                    y,
                                    wachspressKappa,
                                    wachspressA,
                                    wachspressB,
                                    nEdgesOnCellSubset,
                                    vertexIndexSubset):

    numerator = np.zeros((np.shape(x)[0],nEdgesOnCell))
    denominator = np.zeros(np.shape(x)[0])

    for jVertex in range(0,nEdgesOnCell):

        numerator[:,jVertex] = wachspress_numerator(nEdgesOnCell,
                                                    jVertex,
                                                    iVertex,
                                                    x,
                                                    y,
                                                    wachspressKappa,
                                                    wachspressA,
                                                    wachspressB,
                                                    nEdgesOnCellSubset,
                                                    vertexIndexSubset)

        denominator[:] = denominator[:] + numerator[:,jVertex]

    wachpress = np.zeros(np.shape(x)[0])
    wachpress[:] = numerator[:,iVertex] / denominator[:]

    return wachpress

#--------------------------------------------------------
# pwl
#--------------------------------------------------------

def calc_pwl_coefficients(nCells,
                          maxEdges,
                          nEdgesOnCell,
                          xLocal,
                          yLocal,
                          edgesOnCell,
                          dvEdge,
                          areaCell):

    basisSubArea = np.zeros((nCells,maxEdges))
    subBasisGradientU = np.zeros((nCells,maxEdges,3))
    subBasisGradientV = np.zeros((nCells,maxEdges,3))
    subBasisConstant = np.zeros((nCells,maxEdges,3))

    # loop over cells
    for iCell in range(0,nCells):

        alphaPWL = 1.0 / float(nEdgesOnCell[iCell])

        # determine cell centre for piecewise linear basis
        xPWLCentre = 0.0
        yPWLCentre = 0.0

        for iVertexOnCell in range(0, nEdgesOnCell[iCell]):

            xPWLCentre = xPWLCentre + alphaPWL * xLocal[iVertexOnCell,iCell]
            yPWLCentre = yPWLCentre + alphaPWL * yLocal[iVertexOnCell,iCell]

        # calculate the area of the subcells
        basisSubAreaSum = 0.0

        for iSubCell in range(0, nEdgesOnCell[iCell]):

            iEdge = edgesOnCell[iCell,iSubCell]
            iVertexOnCell1 = wrapped_index(iSubCell - 1, nEdgesOnCell[iCell])
            iVertexOnCell2 = iSubCell

            c = dvEdge[iEdge]
            a = math.sqrt(math.pow(xLocal[iVertexOnCell1,iCell] - xPWLCentre,2) + \
                          math.pow(yLocal[iVertexOnCell1,iCell] - yPWLCentre,2))
            b = math.sqrt(math.pow(xLocal[iVertexOnCell2,iCell] - xPWLCentre,2) + \
                          math.pow(yLocal[iVertexOnCell2,iCell] - yPWLCentre,2))

            s = (a + b + c) * 0.5

            # Heron's formula
            basisSubArea[iCell,iSubCell] = math.sqrt(s * (s-a) * (s-b) * (s-c))

            basisSubAreaSum = basisSubAreaSum + basisSubArea[iCell,iSubCell]

        # ensure sum of subareas equals the area of the cell
        basisSubArea[iCell,:] = basisSubArea[iCell,:] * (areaCell[iCell] / basisSubAreaSum)

        # calculate the linear basis on the sub triangle
        for iSubCell in range(0, nEdgesOnCell[iCell]):

            iVertexOnCell1 = wrapped_index(iSubCell - 1, nEdgesOnCell[iCell])
            iVertexOnCell2 = iSubCell

            # set up left hand matrix
            leftMatrix = np.zeros((3,3))

            leftMatrix[0,0] = xLocal[iVertexOnCell1,iCell] - xPWLCentre
            leftMatrix[0,1] = yLocal[iVertexOnCell1,iCell] - yPWLCentre
            leftMatrix[0,2] = 1.0

            leftMatrix[1,0] = xLocal[iVertexOnCell2,iCell] - xPWLCentre
            leftMatrix[1,1] = yLocal[iVertexOnCell2,iCell] - yPWLCentre
            leftMatrix[1,2] = 1.0

            leftMatrix[2,0] = 0.0
            leftMatrix[2,1] = 0.0
            leftMatrix[2,2] = 1.0

            # first basis
            rightHandSide = np.zeros(3)
            rightHandSide[0] = 1.0
            rightHandSide[1] = 0.0
            rightHandSide[2] = 0.0

            lu, piv = lu_factor(leftMatrix)
            solutionVector = lu_solve((lu, piv), rightHandSide)

            subBasisGradientU[iCell,iSubCell,0] = solutionVector[0]
            subBasisGradientV[iCell,iSubCell,0] = solutionVector[1]
            subBasisConstant[iCell,iSubCell,0]  = solutionVector[2]

            # second basis
            rightHandSide[0] = 0.0
            rightHandSide[1] = 1.0
            rightHandSide[2] = 0.0

            lu, piv = lu_factor(leftMatrix)
            solutionVector = lu_solve((lu, piv), rightHandSide)

            subBasisGradientU[iCell,iSubCell,1] = solutionVector[0]
            subBasisGradientV[iCell,iSubCell,1] = solutionVector[1]
            subBasisConstant[iCell,iSubCell,1]  = solutionVector[2]

            # third basis
            subBasisGradientU[iCell,iSubCell,2] = -subBasisGradientU[iCell,iSubCell,0] - subBasisGradientU[iCell,iSubCell,1]
            subBasisGradientV[iCell,iSubCell,2] = -subBasisGradientV[iCell,iSubCell,0] - subBasisGradientV[iCell,iSubCell,1]
            subBasisConstant[iCell,iSubCell,2]  = 1.0 - subBasisConstant[iCell,iSubCell,0] - subBasisConstant[iCell,iSubCell,1]

    return basisSubArea, subBasisGradientU, subBasisGradientV, subBasisConstant

#--------------------------------------------------------

def PointInsideTriangle2(pt,tri):
    '''checks if point pt(2) is inside triangle tri(3x2). @Developer'''
    a = 1/(-tri[1,1]*tri[2,0]+tri[0,1]*(-tri[1,0]+tri[2,0])+ \
        tri[0,0]*(tri[1,1]-tri[2,1])+tri[1,0]*tri[2,1])
    s = a*(tri[2,0]*tri[0,1]-tri[0,0]*tri[2,1]+(tri[2,1]-tri[0,1])*pt[0]+ \
        (tri[0,0]-tri[2,0])*pt[1])
    if s<0: return False
    else: t = a*(tri[0,0]*tri[1,1]-tri[1,0]*tri[0,1]+(tri[0,1]-tri[1,1])*pt[0]+ \
              (tri[1,0]-tri[0,0])*pt[1])
    return ((t>0) and (1-s-t>0))

#--------------------------------------------------------

def pwl_basis_function(nEdgesOnCell,
                       xLocal,
                       yLocal,
                       iBasisVertex,
                       subBasisGradientU,
                       subBasisGradientV,
                       subBasisConstant,
                       x,
                       y):

    pwl = 0.0

    alphaPWL = 1.0 / float(nEdgesOnCell)

    xC = np.mean(xLocal[0:nEdgesOnCell])
    yC = np.mean(yLocal[0:nEdgesOnCell])

    for iVertexOnCell in range(0,nEdgesOnCell):
        iVertexOnCell1 = wrapped_index(iVertexOnCell - 1, nEdgesOnCell)
        iVertexOnCell2 = iVertexOnCell

        if (PointInsideTriangle2(np.array([x,y]),
                                 np.array([[xC,yC],
                                           [xLocal[iVertexOnCell1],yLocal[iVertexOnCell1]],
                                           [xLocal[iVertexOnCell2],yLocal[iVertexOnCell2]]]))):

            pwl += (subBasisConstant[iVertexOnCell,2] + x * subBasisGradientU[iVertexOnCell,2] + y * subBasisGradientV[iVertexOnCell,2]) * alphaPWL

            if (iVertexOnCell == wrapped_index(iBasisVertex + 1, nEdgesOnCell)):

                pwl += subBasisConstant[iVertexOnCell,0] + x * subBasisGradientU[iVertexOnCell,0] + y * subBasisGradientV[iVertexOnCell,0]

            elif (iVertexOnCell == iBasisVertex):

                pwl += subBasisConstant[iVertexOnCell,1] + x * subBasisGradientU[iVertexOnCell,1] + y * subBasisGradientV[iVertexOnCell,1]

    return pwl

#--------------------------------------------------------

def pwl_basis_function_array(nEdgesOnCell,
                             iBasisVertex,
                             iSubCell,
                             subBasisGradientU,
                             subBasisGradientV,
                             subBasisConstant,
                             x,
                             y):

    pwl = np.zeros(np.shape(x)[0])

    alphaPWL = 1.0 / float(nEdgesOnCell)

    pwl[:] += (subBasisConstant[iSubCell,2] + x[:] * subBasisGradientU[iSubCell,2] + y[:] * subBasisGradientV[iSubCell,2]) * alphaPWL

    if (iSubCell == wrapped_index(iBasisVertex + 1, nEdgesOnCell)):

        pwl[:] += subBasisConstant[iSubCell,0] + x[:] * subBasisGradientU[iSubCell,0] + y[:] * subBasisGradientV[iSubCell,0]

    elif (iSubCell == iBasisVertex):

        pwl[:] += subBasisConstant[iSubCell,1] + x[:] * subBasisGradientU[iSubCell,1] + y[:] * subBasisGradientV[iSubCell,1]

    return pwl

#--------------------------------------------------------
# Testing
#--------------------------------------------------------

def plot_basis_function():

    filenameGrid = "grid_hex_0082x0094.nc"

    fileGrid = Dataset(filenameGrid,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    edgesOnCell = fileGrid.variables["edgesOnCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    dvEdge = fileGrid.variables["dvEdge"][:]
    areaCell = fileGrid.variables["areaCell"][:]

    fileGrid.close()

    verticesOnCell[:] -= 1
    edgesOnCell[:] -= 1

    n = 50

    iCell = 123
    iVertexOnCellTest = 5

    xLocal, yLocal = local_coords(nCells,
                                  maxEdges,
                                  nEdgesOnCell,
                                  verticesOnCell,
                                  xCell,
                                  yCell,
                                  xVertex,
                                  yVertex)

    nEdgesOnCellSubset, vertexIndexSubset = wachspress_indexes(nEdgesOnCell[iCell])

    wachspressKappa, wachspressA, wachspressB = calc_wachspress_coefficients(nCells,
                                                                             maxEdges,
                                                                             nEdgesOnCell,
                                                                             xLocal,
                                                                             yLocal)

    basisSubArea, subBasisGradientU, subBasisGradientV, subBasisConstant = calc_pwl_coefficients(nCells,
                                                                                                 maxEdges,
                                                                                                 nEdgesOnCell,
                                                                                                 xLocal,
                                                                                                 yLocal,
                                                                                                 edgesOnCell,
                                                                                                 dvEdge,
                                                                                                 areaCell)



    xvertices = []
    yvertices = []
    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
        iVertex = verticesOnCell[iCell,iVertexOnCell]
        xvertices.append(xVertex[iVertex])
        yvertices.append(yVertex[iVertex])
    iVertex = verticesOnCell[iCell,0]
    xvertices.append(xVertex[iVertex])
    yvertices.append(yVertex[iVertex])

    patches = []
    colours1 = []
    colours2 = []
    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

        iVertexOnCell1 = iVertexOnCell
        iVertexOnCell2 = wrapped_index(iVertexOnCell+1,nEdgesOnCell[iCell])

        iVertex1 = verticesOnCell[iCell,iVertexOnCell1]
        iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

        r1x = xVertex[iVertex1] - xCell[iCell]
        r1y = yVertex[iVertex1] - yCell[iCell]

        r2x = xVertex[iVertex2] - xCell[iCell]
        r2y = yVertex[iVertex2] - yCell[iCell]

        dr1x = r1x / float(n)
        dr1y = r1y / float(n)

        dr2x = r2x / float(n)
        dr2y = r2y / float(n)

        for i in range(0,n):
            for j in range(0,n):

                if (i+j < n):

                    r1xp = xCell[iCell] + float(i) * dr1x + float(j) * dr2x
                    r1yp = yCell[iCell] + float(i) * dr1y + float(j) * dr2y

                    x1 = r1xp + (2.0/6.0) * (dr1x+dr2x)
                    y1 = r1yp + (2.0/6.0) * (dr1y+dr2y)

                    x1 -= xCell[iCell]
                    y1 -= yCell[iCell]

                    wachspressBasis = wachspress_basis_function(nEdgesOnCell[iCell],
                                                                iVertexOnCellTest,
                                                                x1,
                                                                y1,
                                                                wachspressKappa[:,:,iCell],
                                                                wachspressA[:,iCell],
                                                                wachspressB[:,iCell],
                                                                nEdgesOnCellSubset,
                                                                vertexIndexSubset)
                    colours1.append(wachspressBasis)

                    pwl = pwl_basis_function(nEdgesOnCell[iCell],
                                             xLocal[:,iCell],
                                             yLocal[:,iCell],
                                             iVertexOnCellTest,
                                             subBasisGradientU[iCell,:,:],
                                             subBasisGradientV[iCell,:,:],
                                             subBasisConstant[iCell,:,:],
                                             x1,
                                             y1)
                    colours2.append(pwl)

                    patches.append(Polygon([[r1xp,r1yp],
                                            [r1xp+dr1x,r1yp+dr1y],
                                            [r1xp+dr2x,r1yp+dr2y]]))

                    if (i+j < n-1):

                        x2 = r1xp + (2.0/3.0) * (dr1x+dr2x)
                        y2 = r1yp + (2.0/3.0) * (dr1y+dr2y)

                        x2 -= xCell[iCell]
                        y2 -= yCell[iCell]

                        wachspressBasis = wachspress_basis_function(nEdgesOnCell[iCell],
                                                                    iVertexOnCellTest,
                                                                    x2,
                                                                    y2,
                                                                    wachspressKappa[:,:,iCell],
                                                                    wachspressA[:,iCell],
                                                                    wachspressB[:,iCell],
                                                                    nEdgesOnCellSubset,
                                                                    vertexIndexSubset)
                        colours1.append(wachspressBasis)

                        pwl = pwl_basis_function(nEdgesOnCell[iCell],
                                                 xLocal[:,iCell],
                                                 yLocal[:,iCell],
                                                 iVertexOnCellTest,
                                                 subBasisGradientU[iCell,:,:],
                                                 subBasisGradientV[iCell,:,:],
                                                 subBasisConstant[iCell,:,:],
                                                 x2,
                                                 y2)
                        colours2.append(pwl)

                        patches.append(Polygon([[r1xp+dr1x,r1yp+dr1y],
                                                [r1xp+dr2x,r1yp+dr2y],
                                                [r1xp+dr1x+dr2x,r1yp+dr1y+dr2y]]))


    # plot
    fig, axes = plt.subplots(2,1,figsize=(10,20))


    # wachspress
    pc1 = PatchCollection(patches, cmap=plt.get_cmap("jet"))
    pc1.set_array(np.array(colours1))
    axes[0].add_collection(pc1)

    axes[0].plot(xvertices,yvertices,c="k")

    iVertex = verticesOnCell[iCell,iVertexOnCellTest]
    axes[0].scatter(xVertex[iVertex],yVertex[iVertex])

    axes[0].set_aspect("equal")
    axes[0].set_title("(a) Wachspress")

    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(pc1, cax=cax)


    # pwl
    pc2 = PatchCollection(patches, cmap=plt.get_cmap("jet"))
    pc2.set_array(np.array(colours2))
    axes[1].add_collection(pc2)

    axes[1].plot(xvertices,yvertices,c="k")

    iVertex = verticesOnCell[iCell,iVertexOnCellTest]
    axes[1].scatter(xVertex[iVertex],yVertex[iVertex])

    axes[1].set_aspect("equal")
    axes[1].set_title("(b) PWL")

    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(pc2, cax=cax)


    plt.tight_layout()
    plt.savefig("basis.png")

#--------------------------------------------------------

if __name__ == "__main__":

    plot_basis_function()
