import sys
from netCDF4 import Dataset
from math import sqrt, pow, acos, pi, fabs, radians
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from variational_basis_functions import \
    local_coords, \
    wrapped_index, \
    wachspress_indexes, \
    calc_wachspress_coefficients, \
    wachspress_basis_function, \
    wachspress_basis_function_array, \
    calc_pwl_coefficients, \
    pwl_basis_function, \
    pwl_basis_function_array, \
    local_to_global_coordinates, \
    latlon_from_xyz, \
    grid_rotation_forward, \
    grid_rotation_backward, \
    local_eastern_and_northern_unit_vectors
sys.path.append("./strain")
from create_ic import velocities_strains_analytical
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D

#--------------------------------------------------------
# integration
#--------------------------------------------------------

def integration_weights_triangle(integrationOrder):

    # D. A. Dunavant, High degree efficient symmetrical Gaussian quadrature rules for the triangle,
    # Int. J. Num. Meth. Engng, 21(1985):1129-1148.

    if (integrationOrder == 1):

        nIntegrationPoints = 1

        u = [0.33333333333333]

        v = [0.33333333333333]

        weights = [1.00000000000000]

    elif (integrationOrder == 2):

        nIntegrationPoints = 3

        u = [0.16666666666667, 0.16666666666667, 0.66666666666667]

        v = [0.16666666666667, 0.66666666666667, 0.16666666666667]

        weights = [0.33333333333333, 0.33333333333333, 0.33333333333333]

    elif (integrationOrder == 3):

        nIntegrationPoints = 4

        u = [0.33333333333333, 0.20000000000000, 0.20000000000000, 0.60000000000000]

        v = [0.33333333333333, 0.20000000000000, 0.60000000000000, 0.20000000000000]

        weights = [-0.56250000000000, 0.52083333333333, 0.52083333333333, 0.52083333333333]

    elif (integrationOrder == 4):

        nIntegrationPoints = 6

        u = [0.44594849091597, 0.44594849091597, 0.10810301816807, 0.09157621350977, \
             0.09157621350977, 0.81684757298046]

        v = [0.44594849091597, 0.10810301816807, 0.44594849091597, 0.09157621350977, \
             0.81684757298046, 0.09157621350977]

        weights = [0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, \
                   0.10995174365532, 0.10995174365532]

    elif (integrationOrder == 5):

        nIntegrationPoints = 7

        u = [0.33333333333333, 0.47014206410511, 0.47014206410511, 0.05971587178977, \
             0.10128650732346, 0.10128650732346, 0.79742698535309]

        v = [0.33333333333333, 0.47014206410511, 0.05971587178977, 0.47014206410511, \
             0.10128650732346, 0.79742698535309, 0.10128650732346]

        weights = [0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, \
                   0.12593918054483, 0.12593918054483, 0.12593918054483]

    elif (integrationOrder == 6):

        nIntegrationPoints = 12

        u = [0.24928674517091, 0.24928674517091, 0.50142650965818, 0.06308901449150, \
             0.06308901449150, 0.87382197101700, 0.31035245103378, 0.63650249912140, \
             0.05314504984482, 0.63650249912140, 0.31035245103378, 0.05314504984482]

        v = [0.24928674517091, 0.50142650965818, 0.24928674517091, 0.06308901449150, \
             0.87382197101700, 0.06308901449150, 0.63650249912140, 0.05314504984482, \
             0.31035245103378, 0.31035245103378, 0.05314504984482, 0.63650249912140]

        weights = [0.11678627572638, 0.11678627572638, 0.11678627572638, 0.05084490637021, \
                   0.05084490637021, 0.05084490637021, 0.08285107561837, 0.08285107561837, \
                   0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837]

    elif (integrationOrder == 7):

        nIntegrationPoints = 13

        u = [0.33333333333333, 0.26034596607904, 0.26034596607904, 0.47930806784192, \
             0.06513010290222, 0.06513010290222, 0.86973979419557, 0.31286549600487, \
             0.63844418856981, 0.04869031542532, 0.63844418856981, 0.31286549600487, \
             0.04869031542532]

        v = [0.33333333333333, 0.26034596607904, 0.47930806784192, 0.26034596607904, \
             0.06513010290222, 0.86973979419557, 0.06513010290222, 0.63844418856981, \
             0.04869031542532, 0.31286549600487, 0.31286549600487, 0.04869031542532, \
             0.63844418856981]

        weights = [-0.14957004446768, 0.17561525743321, 0.17561525743321, 0.17561525743321, \
                   0.05334723560884, 0.05334723560884, 0.05334723560884, 0.07711376089026, \
                   0.07711376089026, 0.07711376089026, 0.07711376089026, 0.07711376089026, \
                   0.07711376089026]

    elif (integrationOrder == 8):

        nIntegrationPoints = 16

        u = [0.33333333333333, 0.45929258829272, 0.45929258829272, 0.08141482341455, \
             0.17056930775176, 0.17056930775176, 0.65886138449648, 0.05054722831703, \
             0.05054722831703, 0.89890554336594, 0.26311282963464, 0.72849239295540, \
             0.00839477740996, 0.72849239295540, 0.26311282963464, 0.00839477740996]

        v = [0.33333333333333, 0.45929258829272, 0.08141482341455, 0.45929258829272, \
             0.17056930775176, 0.65886138449648, 0.17056930775176, 0.05054722831703, \
             0.89890554336594, 0.05054722831703, 0.72849239295540, 0.00839477740996, \
             0.26311282963464, 0.26311282963464, 0.00839477740996, 0.72849239295540]

        weights = [0.14431560767779, 0.09509163426728, 0.09509163426728, 0.09509163426728, \
                   0.10321737053472, 0.10321737053472, 0.10321737053472, 0.03245849762320, \
                   0.03245849762320, 0.03245849762320, 0.02723031417443, 0.02723031417443, \
                   0.02723031417443, 0.02723031417443, 0.02723031417443, 0.02723031417443]

    elif (integrationOrder == 9):

        nIntegrationPoints = 19

        u = [0.333333333333333, 0.020634961602525, 0.489682519198738, 0.489682519198738, \
             0.125820817014127, 0.437089591492937, 0.437089591492937, 0.623592928761935, \
             0.188203535619033, 0.188203535619033, 0.910540973211095, 0.044729513394453, \
             0.044729513394453, 0.036838412054736, 0.221962989160766, 0.036838412054736, \
             0.741198598784498, 0.221962989160766, 0.741198598784498]

        v = [0.333333333333333, 0.489682519198738, 0.020634961602525, 0.489682519198738, \
             0.437089591492937, 0.125820817014127, 0.437089591492937, 0.188203535619033, \
             0.623592928761935, 0.188203535619033, 0.044729513394453, 0.910540973211095, \
             0.044729513394453, 0.221962989160766, 0.036838412054736, 0.741198598784498, \
             0.036838412054736, 0.741198598784498, 0.221962989160766]

        weights = [0.097135796282799, 0.031334700227139, 0.031334700227139, 0.031334700227139, \
                   0.077827541004774, 0.077827541004774, 0.077827541004774, 0.079647738927210, \
                   0.079647738927210, 0.079647738927210, 0.025577675658698, 0.025577675658698, \
                   0.025577675658698, 0.043283539377289, 0.043283539377289, 0.043283539377289, \
                   0.043283539377289, 0.043283539377289, 0.043283539377289]

    elif (integrationOrder == 10):

        nIntegrationPoints = 25

        u = [0.333333333333333, 0.028844733232685, 0.485577633383657, 0.485577633383657, \
             0.781036849029926, 0.109481575485037, 0.109481575485037, 0.141707219414880, \
             0.307939838764121, 0.141707219414880, 0.550352941820999, 0.307939838764121, \
             0.550352941820999, 0.025003534762686, 0.246672560639903, 0.025003534762686, \
             0.728323904597411, 0.246672560639903, 0.728323904597411, 0.009540815400299, \
             0.066803251012200, 0.009540815400299, 0.923655933587500, 0.066803251012200, \
             0.923655933587500]

        v = [0.333333333333333, 0.485577633383657, 0.028844733232685, 0.485577633383657, \
             0.109481575485037, 0.781036849029926, 0.109481575485037, 0.307939838764121, \
             0.141707219414880, 0.550352941820999, 0.141707219414880, 0.550352941820999, \
             0.307939838764121, 0.246672560639903, 0.025003534762686, 0.728323904597411, \
             0.025003534762686, 0.728323904597411, 0.246672560639903, 0.066803251012200, \
             0.009540815400299, 0.923655933587500, 0.009540815400299, 0.923655933587500, \
             0.066803251012200]

        weights = [0.090817990382754, 0.036725957756467, 0.036725957756467, 0.036725957756467, \
                   0.045321059435528, 0.045321059435528, 0.045321059435528, 0.072757916845420, \
                   0.072757916845420, 0.072757916845420, 0.072757916845420, 0.072757916845420, \
                   0.072757916845420, 0.028327242531057, 0.028327242531057, 0.028327242531057, \
                   0.028327242531057, 0.028327242531057, 0.028327242531057, 0.009421666963733, \
                   0.009421666963733, 0.009421666963733, 0.009421666963733, 0.009421666963733, \
                   0.009421666963733]

    elif (integrationOrder == 12):

        nIntegrationPoints = 33

        u = [0.023565220452390, 0.488217389773805, 0.488217389773805, 0.120551215411079, \
             0.439724392294460, 0.439724392294460, 0.457579229975768, 0.271210385012116, \
             0.271210385012116, 0.744847708916828, 0.127576145541586, 0.127576145541586, \
             0.957365299093579, 0.021317350453210, 0.021317350453210, 0.115343494534698, \
             0.275713269685514, 0.115343494534698, 0.608943235779788, 0.275713269685514, \
             0.608943235779788, 0.022838332222257, 0.281325580989940, 0.022838332222257, \
             0.695836086787803, 0.281325580989940, 0.695836086787803, 0.025734050548330, \
             0.116251915907597, 0.025734050548330, 0.858014033544073, 0.116251915907597, \
             0.858014033544073]

        v = [0.488217389773805, 0.023565220452390, 0.488217389773805, 0.439724392294460, \
             0.120551215411079, 0.439724392294460, 0.271210385012116, 0.457579229975768, \
             0.271210385012116, 0.127576145541586, 0.744847708916828, 0.127576145541586, \
             0.021317350453210, 0.957365299093579, 0.021317350453210, 0.275713269685514, \
             0.115343494534698, 0.608943235779788, 0.115343494534698, 0.608943235779788, \
             0.275713269685514, 0.281325580989940, 0.022838332222257, 0.695836086787803, \
             0.022838332222257, 0.695836086787803, 0.281325580989940, 0.116251915907597, \
             0.025734050548330, 0.858014033544073, 0.025734050548330, 0.858014033544073, \
             0.116251915907597]

        weights = [0.025731066440455, 0.025731066440455, 0.025731066440455, 0.043692544538038, \
                   0.043692544538038, 0.043692544538038, 0.062858224217885, 0.062858224217885, \
                   0.062858224217885, 0.034796112930709, 0.034796112930709, 0.034796112930709, \
                   0.006166261051559, 0.006166261051559, 0.006166261051559, 0.040371557766381, \
                   0.040371557766381, 0.040371557766381, 0.040371557766381, 0.040371557766381, \
                   0.040371557766381, 0.022356773202303, 0.022356773202303, 0.022356773202303, \
                   0.022356773202303, 0.022356773202303, 0.022356773202303, 0.017316231108659, \
                   0.017316231108659, 0.017316231108659, 0.017316231108659, 0.017316231108659, \
                   0.017316231108659]

    return nIntegrationPoints, np.array(u), np.array(v), np.array(weights)

#--------------------------------------------------------

def integration_weights_gauss_lobatto(n):

    if (n == 2):

        u = [-1.0,
              1.0]
        weights = [1.0,
                   1.0]
        nIntegrationPoints = 2

    elif (n == 3):

        u = [-1.0,
              0.0,
              1.0]
        weights = [1.0/3.0,
                   4.0/3.0,
                   1.0/3.0]
        nIntegrationPoints = 3

    elif (n == 4):

        u = [-1.0,
             -1.0 / sqrt(5.0),
              1.0 / sqrt(5.0),
              1.0]
        weights = [1.0/6.0,
                   5.0/6.0,
                   5.0/6.0,
                   1.0/6.0]
        nIntegrationPoints = 4

    elif (n == 5):
        u = [-1.0,
             -sqrt(3.0/7.0),
              0.0,
              sqrt(3.0/7.0),
              1.0]
        weights = [ 1.0/10.0,
                   49.0/90.0,
                   32.0/45.0,
                   49.0/90.0,
                    1.0/10.0]
        nIntegrationPoints = 5

    elif (n == 6):
        u = [-sqrt( (1.0/3.0) - ( (2.0*sqrt(7.0)) / 21.0) ),
              sqrt( (1.0/3.0) - ( (2.0*sqrt(7.0)) / 21.0) ),
             -sqrt( (1.0/3.0) + ( (2.0*sqrt(7.0)) / 21.0) ),
              sqrt( (1.0/3.0) + ( (2.0*sqrt(7.0)) / 21.0) ),
             -1.0,
              1.0]
        weights = [(14 + sqrt(7.0)) / 30.0,
                   (14 + sqrt(7.0)) / 30.0,
                   (14 - sqrt(7.0)) / 30.0,
                   (14 - sqrt(7.0)) / 30.0,
                   1.0 / 15.0,
                   1.0 / 15.0]
        nIntegrationPoints = 6

    elif (n == 7):
        u = [0.0,
              sqrt((5.0/11.0) - (2.0/11.0)*sqrt(5.0/3.0)),
             -sqrt((5.0/11.0) - (2.0/11.0)*sqrt(5.0/3.0)),
              sqrt((5.0/11.0) + (2.0/11.0)*sqrt(5.0/3.0)),
             -sqrt((5.0/11.0) + (2.0/11.0)*sqrt(5.0/3.0)),
             -1.0,
              1.0]
        weights = [256.0/525.0,
                   (124.0 + 7.0*sqrt(15.0)) / 350.0,
                   (124.0 + 7.0*sqrt(15.0)) / 350.0,
                   (124.0 - 7.0*sqrt(15.0)) / 350.0,
                   (124.0 - 7.0*sqrt(15.0)) / 350.0,
                   1.0/21.0,
                   1.0/21.0]
        nIntegrationPoints = 7

    return nIntegrationPoints, u, weights

#--------------------------------------------------------

def triangle_area(x1, y1,
                  x2, y2,
                  x3, y3):

    a = sqrt(pow((x2 - x1),2) +
             pow((y2 - y1),2))
    b = sqrt(pow((x3 - x2),2) +
             pow((y3 - y2),2))
    c = sqrt(pow((x1 - x3),2) +
             pow((y1 - y3),2))

    s = 0.5 * (a + b + c)

    return sqrt(s * (s-a) * (s-b) * (s-c))

#--------------------------------------------------------

def vec_mag(v):

    return sqrt(pow(v[0],2) +
                pow(v[1],2) +
                pow(v[2],2))

#--------------------------------------------------------

def dot_product(a, b):

    dp = a[0] * b[0] + \
         a[1] * b[1] + \
         a[2] * b[2]

    return dp

#--------------------------------------------------------

def cross_product(a, b):

    cp = np.zeros(3)

    cp[0] = a[1] * b[2] - a[2] * b[1]
    cp[1] = a[2] * b[0] - a[0] * b[2]
    cp[2] = a[0] * b[1] - a[1] * b[0]

    return cp

#--------------------------------------------------------

def spherical_angle(r0, r1, r2):

    v1 = np.zeros(3)
    v2 = np.zeros(3)

    v1[:] = r1[:] - r0[:]
    v2[:] = r2[:] - r0[:]

    n1 = cross_product(v1, r0)
    n2 = cross_product(v2, r0)

    n1[:] /= vec_mag(n1)
    n2[:] /= vec_mag(n2)

    angle = acos(dot_product(n1,n2))

    return angle

#--------------------------------------------------------

def spherical_triangle_area(a, b, c, radius):

    A = spherical_angle(a, b, c)
    B = spherical_angle(b, c, a)
    C = spherical_angle(c, a, b)

    area = pow(radius,2) * (A + B + C - pi)

    return area

#--------------------------------------------------------

def get_triangle_mapping(x1, y1,
                         x2, y2,
                         u1, v1,
                         u2, v2):

    mapping = np.zeros((2,2))

    mapping[0,0] = (u2*y1 - u1*y2) / (x2*y1 - x1*y2)
    mapping[0,1] = (u1*x2 - u2*x1) / (y1*x2 - y2*x1)

    mapping[1,0] = (v2*y1 - v1*y2) / (x2*y1 - x1*y2)
    mapping[1,1] = (v1*x2 - v2*x1) / (y1*x2 - y2*x1)

    jacobian = mapping[0,0] * mapping[1,1] - mapping[0,1] * mapping[1,0]

    return mapping, jacobian

#--------------------------------------------------------

def use_triangle_mapping(u, v, mapping):

    x = mapping[0,0] * u + mapping[0,1] * v
    y = mapping[1,0] * u + mapping[1,1] * v

    return x, y

#--------------------------------------------------------

def use_triangle_mapping_array(u, v, mapping):

    x = np.zeros(np.shape(u)[0])
    y = np.zeros(np.shape(u)[0])

    x[:] = mapping[0,0] * u[:] + mapping[0,1] * v[:]
    y[:] = mapping[1,0] * u[:] + mapping[1,1] * v[:]

    return x, y

#--------------------------------------------------------

def L2_norm_wachspress_integral(nCells,
                                maxEdges,
                                nEdgesOnCell,
                                verticesOnCell,
                                areaCell,
                                latCell,
                                lonCell,
                                xCell,
                                yCell,
                                zCell,
                                xVertex,
                                yVertex,
                                zVertex,
                                strain11,
                                strain22,
                                strain12,
                                latitudeLimit,
                                radius):

    mu = 3
    lu = 5

    mv = 2
    lv = 4

    rotateCartesianGrid = True

    xLocal, yLocal = local_coords(nCells,
                                  maxEdges,
                                  nEdgesOnCell,
                                  verticesOnCell,
                                  xCell,
                                  yCell,
                                  zCell,
                                  xVertex,
                                  yVertex,
                                  zVertex,
                                  rotateCartesianGrid)

    wachspressKappa, wachspressA, wachspressB = calc_wachspress_coefficients(nCells,
                                                                             maxEdges,
                                                                             nEdgesOnCell,
                                                                             xLocal,
                                                                             yLocal)

    integrationOrder = 8
    nIntegrationPoints, u, v, weights = integration_weights_triangle(integrationOrder)

    normE11  = 0.0
    denomE11 = 0.0
    normE22  = 0.0
    denomE22 = 0.0
    #normE12  = 0.0
    #denomE12 = 0.0

    #testIntegrationPoints = True

    # testing
    #testAreaIntegration = True
    #if (testAreaIntegration):
    #    areaTest = 0.0
    #    strain11[:] = 1.0
    #    strain22[:] = 1.0
    #    strain12[:] = 1.0

    for iCell in tqdm(range(0,nCells)):

        if (fabs(latCell[iCell]) > np.radians(latitudeLimit)):

            #if (testAreaIntegration):
            #    areaTest += areaCell[iCell]

            nEdgesOnCellSubset, vertexIndexSubset = wachspress_indexes(nEdgesOnCell[iCell])

            xCellr, yCellr, zCellr = grid_rotation_forward(xCell[iCell], yCell[iCell], zCell[iCell], rotateCartesianGrid)

            unitVectorEast, unitVectorNorth = local_eastern_and_northern_unit_vectors(xCellr, yCellr, zCellr)

            #if (testIntegrationPoints):
            #    xp = []
            #    yp = []
            #    zp = []

            for iEdgeOnCell in range(0,nEdgesOnCell[iCell]):

                iEdgeOnCell1 = iEdgeOnCell
                iEdgeOnCell2 = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])

                iVertex1 = verticesOnCell[iCell,iEdgeOnCell1]
                iVertex2 = verticesOnCell[iCell,iEdgeOnCell2]

                areaTriangle = spherical_triangle_area(np.array([xCell[iCell],yCell[iCell],zCell[iCell]]),
                                                       np.array([xVertex[iVertex1],yVertex[iVertex1],zVertex[iVertex1]]),
                                                       np.array([xVertex[iVertex2],yVertex[iVertex2],zVertex[iVertex2]]),
                                                       radius)
                #areaTriangle = triangle_area(0.0, 0.0, \
                #                             xLocal[iEdgeOnCell1,iCell], yLocal[iEdgeOnCell1,iCell], \
                #                             xLocal[iEdgeOnCell2,iCell], yLocal[iEdgeOnCell2,iCell])

                mapping, jacobian = get_triangle_mapping(1.0, 0.0,
                                                         0.0, 1.0,
                                                         xLocal[iEdgeOnCell1,iCell], yLocal[iEdgeOnCell1,iCell],
                                                         xLocal[iEdgeOnCell2,iCell], yLocal[iEdgeOnCell2,iCell])

                # convert from u,v to x,y local
                x2D, y2D = use_triangle_mapping_array(u, v, mapping)

                basisWeights = np.zeros((nEdgesOnCell[iCell],nIntegrationPoints))
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

                    basisWeights[iVertexOnCell,:] = wachspress_basis_function_array(nEdgesOnCell[iCell],
                                                                                    iVertexOnCell,
                                                                                    x2D[:],
                                                                                    y2D[:],
                                                                                    wachspressKappa[:,:,iCell],
                                                                                    wachspressA[:,iCell],
                                                                                    wachspressB[:,iCell],
                                                                                    nEdgesOnCellSubset,
                                                                                    vertexIndexSubset)


                for iWeight in range(0, nIntegrationPoints):

                    # convert x,y to x,y,x
                    x3Dr, y3Dr, z3Dr = local_to_global_coordinates(xCellr, yCellr, zCellr,
                                                                   unitVectorEast,
                                                                   unitVectorNorth,
                                                                   x2D[iWeight], y2D[iWeight])
                    #x3D, y3D, z3D = grid_rotation_backward(x3Dr, y3Dr, z3Dr, rotateCartesianGrid)

                    #if (testIntegrationPoints):
                    #    xp.append(x3D)
                    #    yp.append(y3D)
                    #    zp.append(z3D)

                    # convert x,y,x to lat,lon
                    lat, lon = latlon_from_xyz(x3Dr, y3Dr, z3Dr, radius)

                    # get analytical strains
                    _, _, e11Analytical, e22Analytical, e12Analytical = velocities_strains_analytical(lat, lon, mu, lu, mv, lv)


                    e11 = 0.0
                    e22 = 0.0
                    #e12 = 0.0
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                        #basisWeight = wachspress_basis_function(nEdgesOnCell[iCell],
                        #                                        iVertexOnCell,
                        #                                        x2D,
                        #                                        y2D,
                        #                                        wachspressKappa[:,:,iCell],
                        #                                        wachspressA[:,iCell],
                        #                                        wachspressB[:,iCell],
                        #                                        nEdgesOnCellSubset,
                        #                                        vertexIndexSubset)
                        e11 += strain11[iCell,iVertexOnCell] * basisWeights[iVertexOnCell,iWeight]
                        e22 += strain22[iCell,iVertexOnCell] * basisWeights[iVertexOnCell,iWeight]
                        #e12 += strain12[iCell,iVertexOnCell] * basisWeight


                    #if (not testAreaIntegration):

                    normE11  += weights[iWeight] * areaTriangle * pow(e11 - e11Analytical,2)
                    denomE11 += weights[iWeight] * areaTriangle * pow(e11Analytical,2)

                    normE22  += weights[iWeight] * areaTriangle * pow(e22 - e22Analytical,2)
                    denomE22 += weights[iWeight] * areaTriangle * pow(e22Analytical,2)

                    #normE12  += weights[iWeight] * areaTriangle * pow(e12 - e12Analytical,2)
                    #denomE12 += weights[iWeight] * areaTriangle * pow(e12Analytical,2)

                    #else:

                    #    normE11  += weights[iWeight] * areaTriangle * e11
                    #    normE22  += weights[iWeight] * areaTriangle * e22
                    #    #normE12  += weights[iWeight] * areaTriangle * pow(e12 - e12Analytical,2)

            #if (testIntegrationPoints):

            #    fig = plt.figure()
            #    axis = fig.add_subplot(111, projection='3d')

            #    axis.scatter(xCell[iCell], yCell[iCell], zCell[iCell], c="red")

            #    xv = []
            #    yv = []
            #    zv = []
            #    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            #        iVertex = verticesOnCell[iCell,iVertexOnCell]
            #        axis.plot([xCell[iCell],xVertex[iVertex]],
            #                  [yCell[iCell],yVertex[iVertex]],
            #                  [zCell[iCell],zVertex[iVertex]], c="g")
            #        xv.append(xVertex[iVertex])
            #        yv.append(yVertex[iVertex])
            #        zv.append(zVertex[iVertex])
            #    iVertex = verticesOnCell[iCell,0]
            #    xv.append(xVertex[iVertex])
            #    yv.append(yVertex[iVertex])
            #    zv.append(zVertex[iVertex])
            #    axis.plot(xv, yv, zv)

            #    for i in range(0,len(xp)):
            #        mag = sqrt(pow(xp[i],2) +
            #                   pow(yp[i],2) +
            #                   pow(zp[i],2))
            #        xp[i] *= radius / mag
            #        yp[i] *= radius / mag
            #        zp[i] *= radius / mag

            #    axis.scatter(xp, yp, zp, c="k")

            #    plt.show()

            #    exit(0)


    #if (testAreaIntegration):
    #    print("Test area:     ", areaTest)
    #    print("Integral area: ", normE11, normE22)
    #    exit(0)

    normE11 = sqrt(normE11 / denomE11)
    normE22 = sqrt(normE22 / denomE22)
    #normE12 = sqrt(normE12 / denomE12)

    return normE11, normE22

#--------------------------------------------------------

def get_norm_wachspress_integral(filenameGrid, filenameIn, strain11Name, strain22Name, strain12Name, latitudeLimit):

    fileGrid = Dataset(filenameGrid,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    areaCell = fileGrid.variables["areaCell"][:]
    latCell = fileGrid.variables["latCell"][:]
    lonCell = fileGrid.variables["lonCell"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    verticesOnCell[:] -= 1

    fileIn = Dataset(filenameIn,"r")

    radius = fileIn.sphere_radius

    strain11 = fileIn.variables[strain11Name][0,:,:]
    strain22 = fileIn.variables[strain22Name][0,:,:]
    strain12 = fileIn.variables[strain12Name][0,:,:]

    fileIn.close()

    normE11, normE22 = L2_norm_wachspress_integral(nCells,
                                                   maxEdges,
                                                   nEdgesOnCell,
                                                   verticesOnCell,
                                                   areaCell,
                                                   latCell,
                                                   lonCell,
                                                   xCell,
                                                   yCell,
                                                   zCell,
                                                   xVertex,
                                                   yVertex,
                                                   zVertex,
                                                   strain11,
                                                   strain22,
                                                   strain12,
                                                   latitudeLimit,
                                                   radius)

    return normE11, normE22

#--------------------------------------------------------

def L2_norm_pwl_integral(nCells,
                         maxEdges,
                         nEdgesOnCell,
                         verticesOnCell,
                         edgesOnCell,
                         dvEdge,
                         areaCell,
                         latCell,
                         xCell,
                         yCell,
                         zCell,
                         xVertex,
                         yVertex,
                         zVertex,
                         strain11,
                         strain22,
                         strain12,
                         latitudeLimit,
                         radius):

    mu = 3
    lu = 5

    mv = 2
    lv = 4

    rotateCartesianGrid = True

    xLocal, yLocal = local_coords(nCells,
                                  maxEdges,
                                  nEdgesOnCell,
                                  verticesOnCell,
                                  xCell,
                                  yCell,
                                  zCell,
                                  xVertex,
                                  yVertex,
                                  zVertex,
                                  rotateCartesianGrid)

    basisSubArea, subBasisGradientU, subBasisGradientV, subBasisConstant = calc_pwl_coefficients(nCells,
                                                                                                 maxEdges,
                                                                                                 nEdgesOnCell,
                                                                                                 xLocal,
                                                                                                 yLocal,
                                                                                                 edgesOnCell,
                                                                                                 dvEdge,
                                                                                                 areaCell)

    integrationOrder = 8
    nIntegrationPoints, u, v, weights = integration_weights_triangle(integrationOrder)

    normE11  = 0.0
    denomE11 = 0.0
    normE22  = 0.0
    denomE22 = 0.0
    #normE12  = 0.0
    #denomE12 = 0.0

    #testIntegrationPoints = True

    # testing
    #testAreaIntegration = True
    #if (testAreaIntegration):
    #    areaTest = 0.0
    #    strain11[:] = 1.0
    #    strain22[:] = 1.0
    #    strain12[:] = 1.0

    for iCell in tqdm(range(0,nCells)):

        if (fabs(latCell[iCell]) > np.radians(latitudeLimit)):

            #if (testAreaIntegration):
            #    areaTest += areaCell[iCell]

            xC = 0.0
            yC = 0.0
            xC3D = 0.0
            yC3D = 0.0
            zC3D = 0.0
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                xC += xLocal[iVertexOnCell,iCell]
                yC += yLocal[iVertexOnCell,iCell]
                xC3D += xVertex[iVertex]
                yC3D += yVertex[iVertex]
                zC3D += zVertex[iVertex]
            xC /= float(nEdgesOnCell[iCell])
            yC /= float(nEdgesOnCell[iCell])
            mag = sqrt(pow(xC3D,2) + \
                       pow(yC3D,2) + \
                       pow(zC3D,2))
            xC3D *= (radius / mag)
            yC3D *= (radius / mag)
            zC3D *= (radius / mag)



            nEdgesOnCellSubset, vertexIndexSubset = wachspress_indexes(nEdgesOnCell[iCell])

            xCellr, yCellr, zCellr = grid_rotation_forward(xCell[iCell], yCell[iCell], zCell[iCell], rotateCartesianGrid)

            unitVectorEast, unitVectorNorth = local_eastern_and_northern_unit_vectors(xCellr, yCellr, zCellr)

            #if (testIntegrationPoints):
            #    xp = []
            #    yp = []
            #    zp = []

            for iEdgeOnCell in range(0,nEdgesOnCell[iCell]):

                iEdgeOnCell1 = iEdgeOnCell
                iEdgeOnCell2 = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])

                iVertex1 = verticesOnCell[iCell,iEdgeOnCell1]
                iVertex2 = verticesOnCell[iCell,iEdgeOnCell2]

                areaTriangle = spherical_triangle_area(np.array([xC3D,yC3D,zC3D]),
                                                       np.array([xVertex[iVertex1],yVertex[iVertex1],zVertex[iVertex1]]),
                                                       np.array([xVertex[iVertex2],yVertex[iVertex2],zVertex[iVertex2]]),
                                                       radius)
                #areaTriangle = triangle_area(xC, yC, \
                #                             xLocal[iEdgeOnCell1,iCell], yLocal[iEdgeOnCell1,iCell], \
                #                             xLocal[iEdgeOnCell2,iCell], yLocal[iEdgeOnCell2,iCell])

                mapping, jacobian = get_triangle_mapping(1.0, 0.0,
                                                         0.0, 1.0,
                                                         xLocal[iEdgeOnCell1,iCell]-xC, yLocal[iEdgeOnCell1,iCell]-yC,
                                                         xLocal[iEdgeOnCell2,iCell]-xC, yLocal[iEdgeOnCell2,iCell]-yC)

                # convert from u,v to x,y local
                x2D, y2D = use_triangle_mapping_array(u, v, mapping) # pwl coords

                basisWeights = np.zeros((nEdgesOnCell[iCell],nIntegrationPoints))
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    iSubCell = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])
                    basisWeights[iVertexOnCell,:] = pwl_basis_function_array(nEdgesOnCell[iCell],
                                                                             iVertexOnCell,
                                                                             iSubCell,
                                                                             subBasisGradientU[iCell,:,:],
                                                                             subBasisGradientV[iCell,:,:],
                                                                             subBasisConstant[iCell,:,:],
                                                                             x2D[:]+xC, # xlocal
                                                                             y2D[:]+yC) # ylocal

                for iWeight in range(0, nIntegrationPoints):

                    # convert x,y to x,y,x
                    x3Dr, y3Dr, z3Dr = local_to_global_coordinates(xCellr, yCellr, zCellr,
                                                                   unitVectorEast,
                                                                   unitVectorNorth,
                                                                   x2D[iWeight]+xC, y2D[iWeight]+yC) # x/ylocal
                    #x3D, y3D, z3D = grid_rotation_backward(x3Dr, y3Dr, z3Dr, rotateCartesianGrid)

                    #if (testIntegrationPoints):
                    #    xp.append(x3D)
                    #    yp.append(y3D)
                    #    zp.append(z3D)

                    # convert x,y,x to lat,lon
                    lat, lon = latlon_from_xyz(x3Dr, y3Dr, z3Dr, radius)

                    # get analytical strains
                    _, _, e11Analytical, e22Analytical, e12Analytical = velocities_strains_analytical(lat, lon, mu, lu, mv, lv)

                    e11 = 0.0
                    e22 = 0.0
                    #e12 = 0.0
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                        #basisWeight = pwl_basis_function(nEdgesOnCell[iCell],
                        #                                 xLocal[:,iCell],
                        #                                 yLocal[:,iCell],
                        #                                 iVertexOnCell,
                        #                                 subBasisGradientU[iCell,:,:],
                        #                                 subBasisGradientV[iCell,:,:],
                        #                                 subBasisConstant[iCell,:,:],
                        #                                 x2D,
                        #                                 y2D)
                        e11 += strain11[iCell,iVertexOnCell] * basisWeights[iVertexOnCell,iWeight]
                        e22 += strain22[iCell,iVertexOnCell] * basisWeights[iVertexOnCell,iWeight]
                        #e12 += strain12[iCell,iVertexOnCell] * basisWeights[iVertexOnCell,iWeight]

                    #if (not testAreaIntegration):

                    normE11  += weights[iWeight] * areaTriangle * pow(e11 - e11Analytical,2)
                    denomE11 += weights[iWeight] * areaTriangle * pow(e11Analytical,2)

                    normE22  += weights[iWeight] * areaTriangle * pow(e22 - e22Analytical,2)
                    denomE22 += weights[iWeight] * areaTriangle * pow(e22Analytical,2)

                    #normE12  += weights[iWeight] * areaTriangle * pow(e12 - e12Analytical,2)
                    #denomE12 += weights[iWeight] * areaTriangle * pow(e12Analytical,2)

                    #else:

                    #    normE11  += weights[iWeight] * areaTriangle * e11
                    #    normE22  += weights[iWeight] * areaTriangle * e22
                    #    #normE12  += weights[iWeight] * areaTriangle * pow(e12 - e12Analytical,2)

            #if (testIntegrationPoints):

            #    fig = plt.figure()
            #    axis = fig.add_subplot(111, projection='3d')

            #    axis.scatter(xCell[iCell], yCell[iCell], zCell[iCell], c="red")

            #    axis.scatter(xC3D, yC3D, zC3D, c="blue")

            #    xv = []
            #    yv = []
            #    zv = []
            #    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            #        iVertex = verticesOnCell[iCell,iVertexOnCell]
            #        axis.plot([xC3D,xVertex[iVertex]],
            #                  [yC3D,yVertex[iVertex]],
            #                  [zC3D,zVertex[iVertex]], c="g")
            #        xv.append(xVertex[iVertex])
            #        yv.append(yVertex[iVertex])
            #        zv.append(zVertex[iVertex])
            #    iVertex = verticesOnCell[iCell,0]
            #    xv.append(xVertex[iVertex])
            #    yv.append(yVertex[iVertex])
            #    zv.append(zVertex[iVertex])
            #    axis.plot(xv, yv, zv)

            #    for i in range(0,len(xp)):
            #        mag = sqrt(pow(xp[i],2) +
            #                   pow(yp[i],2) +
            #                   pow(zp[i],2))
            #        xp[i] *= radius / mag
            #        yp[i] *= radius / mag
            #        zp[i] *= radius / mag

            #    axis.scatter(xp, yp, zp, c="k")

            #    plt.show()

            #    exit(0)

    #if (testAreaIntegration):
    #    print("Test area:     ", areaTest)
    #    print("Integral area: ", normE11, normE22)
    #    exit(0)

    normE11 = sqrt(normE11 / denomE11)
    normE22 = sqrt(normE22 / denomE22)
    #normE12 = sqrt(normE12 / denomE12)

    return normE11, normE22

#--------------------------------------------------------

def get_norm_pwl_integral(filenameGrid, filenameIn, strain11Name, strain22Name, strain12Name, latitudeLimit):

    fileGrid = Dataset(filenameGrid,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    edgesOnCell = fileGrid.variables["edgesOnCell"][:]
    dvEdge = fileGrid.variables["dvEdge"][:]
    areaCell = fileGrid.variables["areaCell"][:]
    latCell = fileGrid.variables["latCell"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    verticesOnCell[:] -= 1
    edgesOnCell[:] -= 1

    fileIn = Dataset(filenameIn,"r")

    radius = fileIn.sphere_radius

    strain11 = fileIn.variables[strain11Name][0,:,:]
    strain22 = fileIn.variables[strain22Name][0,:,:]
    strain12 = fileIn.variables[strain12Name][0,:,:]

    fileIn.close()

    normE11, normE22 = L2_norm_pwl_integral(nCells,
                                            maxEdges,
                                            nEdgesOnCell,
                                            verticesOnCell,
                                            edgesOnCell,
                                            dvEdge,
                                            areaCell,
                                            latCell,
                                            xCell,
                                            yCell,
                                            zCell,
                                            xVertex,
                                            yVertex,
                                            zVertex,
                                            strain11,
                                            strain22,
                                            strain12,
                                            latitudeLimit,
                                            radius)

    return normE11, normE22

#--------------------------------------------------------

def L2_norm_weak_integral(nVertices,
                          vertexDegree,
                          edgesOnVertex,
                          cellsOnEdge,
                          cellsOnVertex,
                          dcEdge,
                          latVertex,
                          xCell,
                          yCell,
                          zCell,
                          xVertex,
                          yVertex,
                          zVertex,
                          strain11,
                          strain22,
                          strain12,
                          latitudeLimit,
                          radius):

    mu = 3
    lu = 5

    mv = 2
    lv = 4

    rotateCartesianGrid = True

    integrationOrder = 7
    nIntegrationPoints, u, weights = integration_weights_gauss_lobatto(integrationOrder)

    normE11  = 0.0
    denomE11 = 0.0
    normE22  = 0.0
    denomE22 = 0.0
    #normE12  = 0.0
    #denomE12 = 0.0

    #testIntegrationPoints = False

    #testAreaIntegration = False
    #if (testAreaIntegration):
    #    lengthTest = 0.0
    #    strain11[:] = 1.0
    #    strain22[:] = 1.0
    #    strain12[:] = 1.0

    for iVertex in tqdm(range(0,nVertices)):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            #if (testIntegrationPoints):
            #    xp = []
            #    yp = []
            #    zp = []

            for iEdgeOnVertex in range(0,vertexDegree):

                iEdge = edgesOnVertex[iVertex,iEdgeOnVertex]

                #if (testAreaIntegration):
                #    lengthTest += dcEdge[iEdge]

                iCell1 = cellsOnEdge[iEdge,0]
                iCell2 = cellsOnEdge[iEdge,1]

                for iWeight in range(0, nIntegrationPoints):

                    w1 = 0.5 * (u[iWeight] + 1.0)
                    w2 = 1.0 - w1

                    e11 = w1 * strain11[iCell1] + w2 * strain11[iCell2]
                    e22 = w1 * strain22[iCell1] + w2 * strain22[iCell2]
                    e12 = w1 * strain12[iCell1] + w2 * strain12[iCell2]

                    x3D = w1 * xCell[iCell1] + w2 * xCell[iCell2]
                    y3D = w1 * yCell[iCell1] + w2 * yCell[iCell2]
                    z3D = w1 * zCell[iCell1] + w2 * zCell[iCell2]
                    mag = sqrt(pow(x3D,2) + \
                               pow(y3D,2) + \
                               pow(z3D,2))
                    x3D *= (radius / mag)
                    y3D *= (radius / mag)
                    z3D *= (radius / mag)

                    #if (testIntegrationPoints):
                    #    xp.append(x3D)
                    #    yp.append(y3D)
                    #    zp.append(z3D)

                    # convert x,y,x to lat,lon
                    x3Dr, y3Dr, z3Dr = grid_rotation_forward(x3D, y3D, z3D, rotateCartesianGrid)
                    lat, lon = latlon_from_xyz(x3Dr, y3Dr, z3Dr, radius)

                    # get analytical strains
                    _, _, e11Analytical, e22Analytical, e12Analytical = velocities_strains_analytical(lat, lon, mu, lu, mv, lv)

                    #if (not testAreaIntegration):

                    normE11  += 0.5 * weights[iWeight] * dcEdge[iEdge] * pow(e11 - e11Analytical,2)
                    denomE11 += 0.5 * weights[iWeight] * dcEdge[iEdge] * pow(e11Analytical,2)

                    normE22  += 0.5 * weights[iWeight] * dcEdge[iEdge] * pow(e22 - e22Analytical,2)
                    denomE22 += 0.5 * weights[iWeight] * dcEdge[iEdge] * pow(e22Analytical,2)

                    #normE12  += 0.5 * weights[iWeight] * dcEdge[iEdge] * pow(e12 - e12Analytical,2)
                    #denomE12 += 0.5 * weights[iWeight] * dcEdge[iEdge] * pow(e12Analytical,2)

                    #else:

                    #    normE11  += 0.5 * weights[iWeight] * dcEdge[iEdge] * e11
                    #    normE22  += 0.5 * weights[iWeight] * dcEdge[iEdge] * e22

            #if (testIntegrationPoints):

            #    fig = plt.figure()
            #    axis = fig.add_subplot(111, projection='3d')

            #    axis.scatter(xVertex[iVertex], yVertex[iVertex], zVertex[iVertex], c="red")

            #    xv = []
            #    yv = []
            #    zv = []
            #    for iCellOnVertex in range(0,vertexDegree):
            #        iCell = cellsOnVertex[iVertex,iCellOnVertex]
            #        axis.plot([xVertex[iVertex],xCell[iCell]],
            #                  [yVertex[iVertex],yCell[iCell]],
            #                  [zVertex[iVertex],zCell[iCell]], c="g")
            #        xv.append(xCell[iCell])
            #        yv.append(yCell[iCell])
            #        zv.append(zCell[iCell])
            #    iCell = cellsOnVertex[iVertex,0]
            #    xv.append(xCell[iCell])
            #    yv.append(yCell[iCell])
            #    zv.append(zCell[iCell])
            #    axis.plot(xv, yv, zv, c="b")

            #    for i in range(0,len(xp)):
            #        mag = sqrt(pow(xp[i],2) +
            #                   pow(yp[i],2) +
            #                   pow(zp[i],2))
            #        xp[i] *= radius / mag
            #        yp[i] *= radius / mag
            #        zp[i] *= radius / mag

            #    axis.scatter(xp, yp, zp, c="k")

            #    plt.show()

            #    exit(0)

    #if (testAreaIntegration):
    #    print("Test length:     ", lengthTest)
    #    print("Integral length: ", normE11, normE22)
    #    exit(0)

    normE11 = sqrt(normE11 / denomE11)
    normE22 = sqrt(normE22 / denomE22)
    #normE12 = sqrt(normE12 / denomE12)

    return normE11, normE22

#--------------------------------------------------------

def get_norm_weak_integral(filenameGrid, filenameIn, strain11Name, strain22Name, strain12Name, latitudeLimit):

    fileGrid = Dataset(filenameGrid,"r")

    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    edgesOnVertex = fileGrid.variables["edgesOnVertex"][:]
    cellsOnEdge = fileGrid.variables["cellsOnEdge"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    dcEdge = fileGrid.variables["dcEdge"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    edgesOnVertex[:] -= 1
    cellsOnEdge[:] -= 1
    cellsOnVertex[:] -= 1

    fileIn = Dataset(filenameIn,"r")

    radius = fileIn.sphere_radius

    strain11 = fileGrid.variables[strain11Name][0,:]
    strain22 = fileGrid.variables[strain22Name][0,:]
    strain12 = fileGrid.variables[strain12Name][0,:]

    fileIn.close()

    normE11, normE22 = L2_norm_weak_integral(nVertices,
                                             vertexDegree,
                                             edgesOnVertex,
                                             cellsOnEdge,
                                             cellsOnVertex,
                                             dcEdge,
                                             latVertex,
                                             xCell,
                                             yCell,
                                             zCell,
                                             xVertex,
                                             yVertex,
                                             zVertex,
                                             strain11,
                                             strain22,
                                             strain12,
                                             latitudeLimit,
                                             radius)

    return normE11, normE22

#--------------------------------------------------------
# strain
#--------------------------------------------------------

def L2_norm_strain_vertex(numerical, analytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit):

    # variational

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            for iCellOnVertex in range(0,vertexDegree):
                iCell2 = cellsOnVertex[iVertex,iCellOnVertex]
                if (iCell == iCell2):
                    area = kiteAreasOnVertex[iVertex,iCellOnVertex]

            if (fabs(latVertex[iVertex]) > radians(latitudeLimit)):

                norm  = norm  + area * pow(numerical[iCell,iVertexOnCell] - analytical[iCell,iVertexOnCell],2)

                denom = denom + area * pow(analytical[iCell,iVertexOnCell],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_vertex(numerical, analytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit):

    # variational
    norm = 0.0

    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]

            if (fabs(latVertex[iVertex]) > radians(latitudeLimit)):
                norm = max(norm, fabs(numerical[iCell,iVertexOnCell] - analytical[iCell,iVertexOnCell]))

    return norm

#--------------------------------------------------------

def get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # variational

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    kiteAreasOnVertex = fileGrid.variables["kiteAreasOnVertex"][:]

    verticesOnCell[:] = verticesOnCell[:] - 1
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    fileGrid.close()

    # IC
    fileIC = Dataset(filenameIC, "r")

    # nVertices
    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    #strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    #strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    strain11VarAnalytical = np.zeros((nCells, maxEdges))
    #strain22VarAnalytical = np.zeros((nCells, maxEdges))
    #strain12VarAnalytical = np.zeros((nCells, maxEdges))
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            strain11VarAnalytical[iCell,iVertexOnCell] = strain11VertexAnalytical[iVertex]
            #strain22VarAnalytical[iCell,iVertexOnCell] = strain22VertexAnalytical[iVertex]
            #strain12VarAnalytical[iCell,iVertexOnCell] = strain12VertexAnalytical[iVertex]

    # results
    fileMPAS = Dataset(filename, "r")

    # nCells, maxEdges
    strain11var = fileMPAS.variables["strain11var"][0,:]
    #strain22var = fileMPAS.variables["strain22var"][0,:]
    #strain12var = fileMPAS.variables["strain12var"][0,:]

    fileMPAS.close()

    if (normType == "l2"):
        normE11 = L2_norm_strain_vertex(strain11var, strain11VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_vertex(strain22var, strain22VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_vertex(strain12var, strain12VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_vertex(strain11var, strain11VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_vertex(strain22var, strain22VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_vertex(strain12var, strain12VarAnalytical, nCells, vertexDegree, nEdgesOnCell, verticesOnCell, cellsOnVertex, kiteAreasOnVertex, latVertex, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_strain_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    # weak

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (fabs(latCell[iCell]) > np.radians(latitudeLimit)):

            norm  = norm  + areaCell[iCell] * pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * pow(analytical[iCell],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_cell(numerical, analytical, nCells, latCell, areaCell, latitudeLimit):

    # weak
    norm  = 0.0
    for iCell in range(0,nCells):

        if (fabs(latCell[iCell]) > np.radians(latitudeLimit)):
            norm = max(norm, fabs(numerical[iCell] - analytical[iCell]))

    return norm

#--------------------------------------------------------

def get_norm_strain_cell(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # weak

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    latCell = fileGrid.variables["latCell"][:]
    areaCell = fileGrid.variables["areaCell"][:]

    fileGrid.close()

    # IC
    fileIC = Dataset(filenameIC, "r")

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    #strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    #strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    #strain22weak = fileMPAS.variables["strain22weak"][0,:]
    #strain12weak = fileMPAS.variables["strain12weak"][0,:]

    if (normType == "l2"):
        normE11 = L2_norm_strain_cell(strain11weak, strain11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_cell(strain22weak, strain22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_cell(strain12weak, strain12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_cell(strain11weak, strain11CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_cell(strain22weak, strain22CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_cell(strain12weak, strain12CellAnalytical, nCells, latCell, areaCell, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    fileMPAS.close()

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_strain_vertex_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # variational averaged

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            norm  = norm  + areaTriangle[iVertex] * pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * pow(analytical[iVertex],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_vertex_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # variational averaged
    norm  = 0.0

    for iVertex in range(0,nVertices):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):
            norm = max(norm, fabs(numerical[iVertex] - analytical[iVertex]))

    return norm

#--------------------------------------------------------

def get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # variational averaged

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    kiteAreasOnVertex = fileGrid.variables["kiteAreasOnVertex"][:]
    areaTriangle = fileGrid.variables["areaTriangle"][:]

    verticesOnCell[:] = verticesOnCell[:] - 1
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    fileGrid.close()

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    #strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    #strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    strain11var = fileMPAS.variables["strain11var"][0,:]
    #strain22var = fileMPAS.variables["strain22var"][0,:]
    #strain12var = fileMPAS.variables["strain12var"][0,:]

    fileMPAS.close()

    strain11varAvg = np.zeros(nVertices)
    #strain22varAvg = np.zeros(nVertices)
    #strain12varAvg = np.zeros(nVertices)
    denom = np.zeros(nVertices)
    for iVertex in range(0,nVertices):
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex2 = verticesOnCell[iCell,iVertexOnCell]
                if (iVertex2 == iVertex):
                    strain11varAvg[iVertex] += strain11var[iCell,iVertexOnCell]
                    #strain22varAvg[iVertex] += strain22var[iCell,iVertexOnCell]
                    #strain12varAvg[iVertex] += strain12var[iCell,iVertexOnCell]
                    denom[iVertex] += 1.0
    for iVertex in range(0,nVertices):
        if (denom[iVertex] > 0.0):
            strain11varAvg[iVertex] /= denom[iVertex]
            #strain22varAvg[iVertex] /= denom[iVertex]
            #strain12varAvg[iVertex] /= denom[iVertex]

    if (normType == "l2"):
        normE11 = L2_norm_strain_vertex_avg(strain11varAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_vertex_avg(strain22varAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_vertex_avg(strain12varAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_vertex_avg(strain11varAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_vertex_avg(strain22varAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_vertex_avg(strain12varAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    return normE11, normE22, normE12

#--------------------------------------------------------

def L2_norm_strain_cell_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # weak averaged

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            norm  = norm  + areaTriangle[iVertex] * pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * pow(analytical[iVertex],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_strain_cell_avg(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    # weak averaged
    norm  = 0.0

    for iVertex in range(0,nVertices):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):
            norm = max(norm, fabs(numerical[iVertex] - analytical[iVertex]))

    return norm

#--------------------------------------------------------

def get_norm_strain_cell_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType):

    # weak averaged

    # grid
    fileGrid = Dataset(filenameGrid, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]
    kiteAreasOnVertex = fileGrid.variables["kiteAreasOnVertex"][:]
    areaCell = fileGrid.variables["areaCell"][:]
    areaTriangle = fileGrid.variables["areaTriangle"][:]

    verticesOnCell[:] = verticesOnCell[:] - 1
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    fileGrid.close()

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    #strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    #strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    #strain22weak = fileMPAS.variables["strain22weak"][0,:]
    #strain12weak = fileMPAS.variables["strain12weak"][0,:]

    fileMPAS.close()

    strain11weakAvg = np.zeros(nVertices)
    #strain22weakAvg = np.zeros(nVertices)
    #strain12weakAvg = np.zeros(nVertices)
    denom = np.zeros(nVertices)
    for iVertex in range(0,nVertices):
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell > 0 and iCell < nCells):
                strain11weakAvg[iVertex] += strain11weak[iCell]
                #strain22weakAvg[iVertex] += strain22weak[iCell]
                #strain12weakAvg[iVertex] += strain12weak[iCell]
                denom[iVertex] += 1.0
    for iVertex in range(0,nVertices):
        if (denom[iVertex] > 0.0):
            strain11weakAvg[iVertex] /= denom[iVertex]
            #strain22weakAvg[iVertex] /= denom[iVertex]
            #strain12weakAvg[iVertex] /= denom[iVertex]

    if (normType == "l2"):
        normE11 = L2_norm_strain_cell_avg(strain11weakAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#L2_norm_strain_cell_avg(strain22weakAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#L2_norm_strain_cell_avg(strain12weakAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    elif (normType == "linf"):
        normE11 = Linf_norm_strain_cell_avg(strain11weakAvg, strain11VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE22 = 0.0#Linf_norm_strain_cell_avg(strain22weakAvg, strain22VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normE12 = 0.0#Linf_norm_strain_cell_avg(strain12weakAvg, strain12VertexAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    return normE11, normE22, normE12

#--------------------------------------------------------
# stress divergence
#--------------------------------------------------------

def L2_norm_stress_divergence(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):

            norm  = norm  + areaTriangle[iVertex] * pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * pow(analytical[iVertex],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def Linf_norm_stress_divergence(numerical, analytical, nVertices, latVertex, areaTriangle, latitudeLimit):

    degreesToRadians = pi / 180.0

    norm  = 0.0

    for iVertex in range(0,nVertices):

        if (fabs(latVertex[iVertex]) > np.radians(latitudeLimit)):
            norm = max(norm, fabs(numerical[iVertex] - analytical[iVertex]))

    return norm

#--------------------------------------------------------

def get_norm_stress_divergence(filenameIC, filename, latitudeLimit, normType):

    fileIC = Dataset(filenameIC, "r")

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    #stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    latVertex = fileMPAS.variables["latVertex"][:]

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    stressDivergenceU = fileMPAS.variables["stressDivergenceU"][0,:]
    #stressDivergenceV = fileMPAS.variables["stressDivergenceV"][0,:]

    if (normType == "l2"):
        normU = L2_norm_stress_divergence(stressDivergenceU, stressDivergenceUAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normV = 0.0#L2_norm_stress_divergence(stressDivergenceV, stressDivergenceVAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    elif (normType == "linf"):
        normU = Linf_norm_stress_divergence(stressDivergenceU, stressDivergenceUAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
        normV = 0.0#Linf_norm_stress_divergence(stressDivergenceV, stressDivergenceVAnalytical, nVertices, latVertex, areaTriangle, latitudeLimit)
    else:
        raise Exception("Unknown normType: " + normType)

    fileMPAS.close()

    return normU, normV

#--------------------------------------------------------

def get_resolution(filename, latitudeLimit):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    degreesToRadians = pi / 180.0

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        if (fabs(latEdge[iEdge]) > latitudeLimit * degreesToRadians):
            resolution = resolution + dcEdge[iEdge]
            denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

#--------------------------------------------------------

def scaling_lines(axis, xMin, xMax, yMin):

    # linear scaling
    scale = yMin / pow(xMin,1)
    scaleMinLin = pow(xMin,1) * scale
    scaleMaxLin = pow(xMax,1) * scale

    # quadratic scaling
    scale = yMin / pow(xMin,2)
    scaleMinQuad = pow(xMin,2) * scale
    scaleMaxQuad = pow(xMax,2) * scale

    axis.loglog([xMin, xMax], [scaleMinLin,  scaleMaxLin],  linestyle=':', color='k')
    axis.loglog([xMin, xMax], [scaleMinQuad, scaleMaxQuad], linestyle=':', color='k')

#--------------------------------------------------------

def strain_scaling(axes, normType, label, xlabel, xMin, xMax, yMin):

    resolutions = [2562,10242,40962,163842]
    #resolutions = [2562]

    methods = ["wachspress","pwl","weak","wachspress_avg","pwl_avg","weak_avg"]
    #methods = ["weak"]

    latitudeLimit = 20.0

    dirname = {"wachspress":"wachspress",
               "pwl":"pwl",
               "weak":"weak",
               "wachspress_avg":"wachspress",
               "pwl_avg":"pwl",
               "weak_avg":"weak"}

    # plot options
    legendLabels = {"wachspress":"Wachs.",
                    "pwl":"PWL",
                    "weak":"FV",
                    "wachspress_avg":"Wachs. Avg.",
                    "pwl_avg":"PWL Avg.",
                    "weak_avg":"FV Avg."}

    lineColors = {"wachspress":"black",
                  "pwl":"grey",
                  "weak":"red",
                  "wachspress_avg":"blue",
                  "pwl_avg":"green",
                  "weak_avg":"darkturquoise"}

    lineStyles = {"wachspress":"solid",
                  "pwl":"solid",
                  "weak":"solid",
                  "wachspress_avg":"solid",
                  "pwl_avg":"solid",
                  "weak_avg":"solid"}

    markers = {"wachspress":"+",
               "pwl":"x",
               "weak":"^",
               "wachspress_avg":"+",
               "pwl_avg":"x",
               "weak_avg":"^"}

    # plot
    # scaling lines
    scaling_lines(axes, xMin, xMax, yMin)

    plotHandles = []

    for method in methods:

        x = []
        y = []

        for resolution in resolutions:

            filenameGrid = "./strain/grid.%i.nc" %(resolution)
            filenameIC = "./strain/ic_%i.nc" %(resolution)
            filename = "./strain/output_%s_%i/output.2000.nc" %(dirname[method],resolution)

            print(filename, filenameIC)

            if (method == "wachspress"):
                if (normType == "l2"):
                    #normE11, normE22, normE12 = get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType)
                    normE11, normE22 = get_norm_wachspress_integral(filenameGrid, filename, "strain11var", "strain22var", "strain12var", latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "pwl"):
                if (normType == "l2"):
                    #normE11, normE22, normE12 = get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType)
                    normE11, normE22 = get_norm_pwl_integral(filenameGrid, filename, "strain11var", "strain22var", "strain12var", latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_strain_vertex(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "weak"):
                if (normType == "l2"):
                    #normE11, normE22, normE12 = get_norm_strain_cell(filenameGrid, filenameIC, filename, latitudeLimit, normType)
                    normE11, normE22 = get_norm_weak_integral(filename, filename, "strain11weak", "strain22weak", "strain12weak", latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_strain_cell(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "wachspress_avg"):
                if (normType == "l2"):
                    #normE11, normE22, normE12 = get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)
                    normE11, normE22 = get_norm_wachspress_integral(filenameGrid, filename, "strain11varAvg", "strain22varAvg", "strain12varAvg", latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "pwl_avg"):
                if (normType == "l2"):
                    #normE11, normE22, normE12 = get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)
                    normE11, normE22 = get_norm_pwl_integral(filenameGrid, filename, "strain11varAvg", "strain22varAvg", "strain12varAvg", latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_strain_vertex_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)
            elif (method == "weak_avg"):
                if (normType == "l2"):
                    #normE11, normE22, normE12 = get_norm_strain_cell_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)
                    normE11, normE22 = get_norm_wachspress_integral(filename, filename, "strain11weakAvg", "strain22weakAvg", "strain12weakAvg", latitudeLimit)
                else:
                    normE11, normE22, normE12 = get_norm_strain_cell_avg(filenameGrid, filenameIC, filename, latitudeLimit, normType)

            x.append(get_resolution(filename, latitudeLimit))
            y.append(normE11)

        print(y)
        plotHandle, = axes.loglog(x,y, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])
        plotHandles.append(plotHandle)

    #axes.legend(frameon=False, loc=2, fontsize=8, handlelength=4)
    legend1 = axes.legend(handles=plotHandles[0:3], frameon=False, loc=2, fontsize=8, handlelength=2)
    legend2 = axes.legend(handles=plotHandles[3:6], frameon=False, loc=4, fontsize=8, handlelength=2)
    axes.add_artist(legend1)

    axes.set_xlabel(xlabel)
    if (normType == "l2"):
        axes.set_ylabel(r"$L_2$ error norm")
    elif (normType == "linf"):
        axes.set_ylabel(r"$L_\infty$ error norm")
    #axes.set_xlim([xMin, xMax])
    ###axes.set_ylim([None, 2e-1])
    axes.set_title(r'%s $\dot{\epsilon}_{11}$' %(label) ,loc="left")
    axes.set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes.set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes.set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes.set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)
    #axes.set_xticklabels(labels=[1e-2,5e-2],minor=False)

#--------------------------------------------------------

def stress_divergence_scaling(axes, normType, label, xlabel, xMin, xMax, yMin):

    resolutions = [2562,10242,40962,163842]
    #resolutions = [2562,10242]

    methods = ["wachspress", "pwl", "weak", "wachspress_alt", "pwl_alt"]

    latitudeLimit = 20.0

    # plot options
    legendLabels = {"wachspress":"Wachs.",
                    "pwl":"PWL",
                    "weak":"FV",
                    "wachspress_alt":"Wachs. alt",
                    "pwl_alt":"PWL alt."}

    lineColors = {"wachspress":"black",
                  "pwl":"grey",
                  "weak":"red",
                  "wachspress_alt":"black",
                  "pwl_alt":"grey"}

    lineStyles = {"wachspress":"solid",
                  "pwl":"solid",
                  "weak":"solid",
                  "wachspress_alt":"dashed",
                  "pwl_alt":"dashed"}

    markers = {"wachspress":"+",
               "pwl":"x",
               "weak":"^",
               "wachspress_alt":"+",
               "pwl_alt":"x"}

    # plot
    # scaling lines
    scaling_lines(axes, xMin, xMax, yMin)

    for method in methods:

        x = []
        y = []

        for resolution in resolutions:

            filename = "./stress_divergence/output_%s_%i/output.2000.nc" %(method,resolution)
            filenameIC = "./stress_divergence/ic_%i.nc" %(resolution)

            print(filename, filenameIC)

            normU, normV = get_norm_stress_divergence(filenameIC, filename, latitudeLimit, normType)

            x.append(get_resolution(filename, latitudeLimit))
            y.append(normU)

        axes.loglog(x,y, marker=markers[method], color=lineColors[method], ls=lineStyles[method], markersize=5.0, label=legendLabels[method])

    axes.legend(frameon=False, loc=4, fontsize=8, handlelength=4)

    axes.set_xlabel(xlabel)
    #axes.set_ylabel(r"$L_2$ error norm")
    axes.set_ylabel(None)
    #axes.set_xlim([8e-3, 8.5e-2])
    axes.set_title(r'%s $(\nabla \cdot \sigma)_u$' %(label), loc="left")
    axes.set_xticks(ticks=[9e-3,2e-2,3e-2,4e-2,6e-2,7e-2,8e-2],minor=True)
    axes.set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
    axes.set_xticks(ticks=[1e-2,5e-2],minor=False)
    axes.set_xticklabels(labels=[r'$10^{-2}$',r'$5\times10^{-2}$'],minor=False)
    #axes.set_xticklabels(labels=[1e-2,5e-2],minor=False)

#-------------------------------------------------------------------------------

def spherical_operators_scaling():

    cm = 1/2.54  # centimeters in inches
    plt.rc('font',family="Times New Roman")
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

    fig, axes = plt.subplots(1, 2, figsize=(15*cm,7*cm))

    strain_scaling(axes[0], "l2", "(a)", "Grid resolution")

    stress_divergence_scaling(axes[1], "l2", "(b)", "Grid resolution")

    #plt.tight_layout()
    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("spherical_operators_scaling.png",dpi=300)
    plt.savefig("spherical_operators_scaling.eps")

#-------------------------------------------------------------------------------

def spherical_operators_scaling_L2_Linf():

    cm = 1/2.54  # centimeters in inches
    plt.rc('font',family="Times New Roman")
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

    fig, axes = plt.subplots(2, 2, figsize=(15*cm,14*cm))

    xMin = 1.5e-2
    xMax = 3e-2
    yMin = 3e-4

    strain_scaling(axes[0,0], "l2", "(a)", None, xMin, xMax, yMin)

    xMin = 1e-2
    xMax = 2e-2
    yMin = 7e-3

    stress_divergence_scaling(axes[0,1], "l2", "(b)", None, xMin, xMax, yMin)

    xMin = 1.5e-2
    xMax = 3e-2
    yMin = 1e-3

    strain_scaling(axes[1,0], "linf", "(c)", "Grid resolution", xMin, xMax, yMin)

    xMin = 1e-2
    xMax = 2e-2
    yMin = 2.2e-1

    stress_divergence_scaling(axes[1,1], "linf", "(d)", "Grid resolution", xMin, xMax, yMin)

    #plt.tight_layout()
    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("spherical_operators_scaling_l2linf.png",dpi=300)
    plt.savefig("spherical_operators_scaling_l2linf.eps")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    spherical_operators_scaling_L2_Linf()
