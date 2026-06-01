import cartopy.crs as ccrs
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from math import radians, degrees
import matplotlib.pyplot as plt
import sys
import matplotlib

#-------------------------------------------------------------------------------

def setup_maps_projection(location):

    # Geographic coordinate system of input data
    src_crs = ccrs.PlateCarree()

    # Stereographic projection centered at desired point
    if (location == "greenland"):
        proj = ccrs.Stereographic(
            central_longitude=-38.45,
            central_latitude=72.58
        )
    elif (location == "antarctica"):
        proj = ccrs.Stereographic(
            central_longitude=0.0,
            central_latitude=-90.0
        )
    else:
        raise Exception("Unknown projection location")

    return src_crs, proj

#-------------------------------------------------------------------------------

def projection(lat,lon,proj,src_crs):

    x, y = proj.transform_point(
        lon,
        lat,
        src_crs
    )

    return x, y

#-------------------------------------------------------------------------------

def projection_list(lat,lon,proj,src_crs):

    xy = proj.transform_points(
        src_crs,
        lon,
        lat
    )

    return xy[:,0], xy[:,1]

#-------------------------------------------------------------------------------

def plot_limits(xMin, xMax, yMin, yMax, enlarge=1.05, square=False):

    dx = xMax-xMin
    dy = yMax-yMin
    dxy = max(dx,dy) * 0.5
    xc = 0.5*(xMin+xMax)
    yc = 0.5*(yMin+yMax)
    if (square):
        dxUse = dxy
        dyUse = dxy
    else:
        dxUse = dx * 0.5
        dyUse = dy * 0.5
    xMinOut = xc - dxUse*enlarge
    xMaxOut = xc + dxUse*enlarge
    yMinOut = yc - dyUse*enlarge
    yMaxOut = yc + dyUse*enlarge

    return xMinOut, xMaxOut, yMinOut, yMaxOut

#-------------------------------------------------------------------------------

def mesh_patches(location,
                 proj,
                 src_crs,
                 nEdges,
                 nCells,
                 nEdgesOnCell,
                 cellsOnEdge,
                 latEdge,
                 latCell,
                 verticesOnEdge,
                 verticesOnCell,
                 latVertex,
                 lonVertex):

    boundaryEdge = np.zeros(nEdges,dtype="i")
    for iEdge in range(0,nEdges):
        if (cellsOnEdge[iEdge,0] == -1 or
            cellsOnEdge[iEdge,1] == -1):
            boundaryEdge[iEdge] = 1

    lineSegments = []
    for iEdge in range(0,nEdges):
        if (boundaryEdge[iEdge] == 1 and
            ((location == "antarctica" and latEdge[iEdge] < radians(-30.0)) or
             (location == "greenland"  and latEdge[iEdge] > radians( 30.0)))):
            iVertex1 = verticesOnEdge[iEdge,0]
            iVertex2 = verticesOnEdge[iEdge,1]
            x1, y1 = projection(degrees(latVertex[iVertex1]),
                                degrees(lonVertex[iVertex1]),
                                proj,src_crs)
            x2, y2 = projection(degrees(latVertex[iVertex2]),
                                degrees(lonVertex[iVertex2]),
                                proj,src_crs)
            lineSegments.append([[x1,y1],
                                 [x2,y2]])

    lc = LineCollection(lineSegments, color="black", linestyle='solid', linewidth=0.2)

    # plot mesh
    patches = []
    for iCell in range(0,nCells):
        if ((location == "antarctica" and latCell[iCell] < radians(-30.0)) or
            (location == "greenland"  and latCell[iCell] > radians( 30.0))):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                x, y = projection(degrees(latVertex[iVertex]),
                                  degrees(lonVertex[iVertex]),
                                  proj,src_crs)
                vertices.append([x,y])
            patches.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))

    pc = PatchCollection(patches, match_original=True)

    return pc, lc

#-------------------------------------------------------------------------------

def plot_cell_field(field,
                    title,
                    cbtitle,
                    filenameOut,
                    location,
                    nEdges,
                    nCells,
                    nEdgesOnCell,
                    cellsOnEdge,
                    latEdge,
                    latCell,
                    verticesOnEdge,
                    verticesOnCell,
                    latVertex,
                    lonVertex):

    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
    })

    src_crs, proj = setup_maps_projection(location)

    # start plot
    fig = plt.figure()
    axis = plt.axes(projection=proj)

    fig, axis, _, _, _, _ = plot_cell_field_axis(fig,
                                                 axis,
                                                 field,
                                                 title,
                                                 cbtitle,
                                                 src_crs,
                                                 proj,
                                                 location,
                                                 nEdges,
                                                 nCells,
                                                 nEdgesOnCell,
                                                 cellsOnEdge,
                                                 latEdge,
                                                 latCell,
                                                 verticesOnEdge,
                                                 verticesOnCell,
                                                 latVertex,
                                                 lonVertex)

    plt.tight_layout()
    plt.savefig(filenameOut,dpi=600)
    plt.close()

#-------------------------------------------------------------------------------

def plot_cell_field_axis(fig,
                         axis,
                         field,
                         title,
                         cbtitle,
                         src_crs,
                         proj,
                         location,
                         nEdges,
                         nCells,
                         nEdgesOnCell,
                         cellsOnEdge,
                         latEdge,
                         latCell,
                         verticesOnEdge,
                         verticesOnCell,
                         latVertex,
                         lonVertex,
                         xMinIn=None,
                         xMaxIn=None,
                         yMinIn=None,
                         yMaxIn=None):

    axis.set_facecolor('grey')

    # plot mesh
    pc, lc = mesh_patches(location,
                          proj,
                          src_crs,
                          nEdges,
                          nCells,
                          nEdgesOnCell,
                          cellsOnEdge,
                          latEdge,
                          latCell,
                          verticesOnEdge,
                          verticesOnCell,
                          latVertex,
                          lonVertex)

    axis.add_collection(pc)
    axis.add_collection(lc)

    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max

    patches = []
    colors = []
    for iCell in range(0,nCells):
        if (((location == "antarctica" and latCell[iCell] < radians(-40.0)) or
             (location == "greenland"  and latCell[iCell] > radians( 40.0))) and
            field[iCell] > 0.0):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                x, y = projection(degrees(latVertex[iVertex]),
                                  degrees(lonVertex[iVertex]),
                                  proj,src_crs)

                xMin = min(xMin,x)
                xMax = max(xMax,x)
                yMin = min(yMin,y)
                yMax = max(yMax,y)

                vertices.append([x,y])
            patches.append(Polygon(vertices, closed=True, edgecolor="grey", facecolor="white", linewidth=0.1))
            colors.append(field[iCell])

    colors = np.array(colors)
    vmax = np.amax(colors)

    xMin, xMax, yMin, yMax = plot_limits(xMin, xMax, yMin, yMax)

    if (xMinIn is not None and
        xMaxIn is not None and
        yMinIn is not None and
        yMaxIn is not None):
        xMinUse = xMinIn
        xMaxUse = xMaxIn
        yMinUse = yMinIn
        yMaxUse = yMaxIn
    else:
        xMinUse = xMin
        xMaxUse = xMax
        yMinUse = yMin
        yMaxUse = yMax

    cmap = matplotlib.colormaps.get_cmap("viridis")
    pc = PatchCollection(patches, match_original=True, cmap=cmap, norm=matplotlib.colors.LogNorm())
    pc.set_array(colors)

    axis.add_collection(pc)

    gl = axis.gridlines(linewidth=0.5,linestyle="dashed",draw_labels=False)

    #axis.autoscale_view()
    axis.set_xlim((xMinUse,xMaxUse))
    axis.set_ylim((yMinUse,yMaxUse))

    axis.set_aspect("equal")
    axis.set_title(title)
    fig.colorbar(pc,label=cbtitle)

    return fig, axis, xMin, xMax, yMin, yMax

#-------------------------------------------------------------------------------

def data_range(field,
               src_crs,
               proj,
               location,
               nCells,
               nEdgesOnCell,
               verticesOnCell,
               latCell,
               latVertex,
               lonVertex):

    xMin =  sys.float_info.max
    xMax = -sys.float_info.max
    yMin =  sys.float_info.max
    yMax = -sys.float_info.max

    for iCell in range(0,nCells):
        if (((location == "antarctica" and latCell[iCell] < radians(-40.0)) or
             (location == "greenland"  and latCell[iCell] > radians( 40.0))) and
            field[iCell] > 0.0):
            vertices = []
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                x, y = projection(degrees(latVertex[iVertex]),
                                  degrees(lonVertex[iVertex]),
                                  proj,src_crs)

                xMin = min(xMin,x)
                xMax = max(xMax,x)
                yMin = min(yMin,y)
                yMax = max(yMax,y)

    return xMin, xMax, yMin, yMax

#-------------------------------------------------------------------------------
