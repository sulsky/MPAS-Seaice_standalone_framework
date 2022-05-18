from netCDF4 import Dataset
from math import sqrt, pow, pi
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
    pwl_basis_function_array
from create_ics import velocities_strains_stress_divergences
from tqdm import tqdm

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

def get_real_coords_from_barycentric_coords_coefficients(x1, y1, x2, y2, x3, y3):

    Cx = x1
    Cy = y1

    T11 = x2 - x1
    T21 = y2 - y1

    T12 = x3 - x1
    T22 = y3 - y1

    return T11, T12, T21, T22, Cx, Cy

#--------------------------------------------------------

def get_real_coords_from_barycentric_coords(u, v, T11, T12, T21, T22, Cx, Cy):

    x = T11 * u + T12 * v + Cx

    y = T21 * u + T22 * v + Cy

    return x, y

#--------------------------------------------------------

def get_real_coords_from_barycentric_coords_array(u, v, T11, T12, T21, T22, Cx, Cy):

    x = np.zeros(u.shape[0])
    y = np.zeros(u.shape[0])

    x[:] = T11 * u[:] + T12 * v[:] + Cx

    y[:] = T21 * u[:] + T22 * v[:] + Cy

    return x, y

#--------------------------------------------------------

def get_analytical_strain_triangle(u, v, T11, T12, T21, T22, Cx, Cy, xMin, yMin):

    x, y = get_real_coords_from_barycentric_coords(u, v, T11, T12, T21, T22, Cx, Cy)

    uVel, vVel, e11, e22, e12, divu, divv = velocities_strains_stress_divergences(x-xMin, y-yMin, "sinusoid")

    return e11, e22, e12, x, y

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

def L2_norm_wachspress_integral(nCells,
                                maxEdges,
                                nEdgesOnCell,
                                verticesOnCell,
                                areaCell,
                                xCell,
                                yCell,
                                xVertex,
                                yVertex,
                                useCell,
                                strain11,
                                strain22,
                                strain12):

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    xLocal, yLocal = local_coords(nCells,
                                  maxEdges,
                                  nEdgesOnCell,
                                  verticesOnCell,
                                  xCell,
                                  yCell,
                                  xVertex,
                                  yVertex)

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

    # testing
    #testIntegrationPoints = True
    #testAreaIntegration = False
    #if (testAreaIntegration):
    #    areaTest = 0.0
    #    strain11[:] = 1.0
    #    strain22[:] = 1.0
    #    strain12[:] = 1.0

    for iCell in tqdm(range(0,nCells)):

        if (useCell[iCell]):

            #if (testAreaIntegration):
            #    areaTest += areaCell[iCell]

            nEdgesOnCellSubset, vertexIndexSubset = wachspress_indexes(nEdgesOnCell[iCell])

            #if (testIntegrationPoints):
            #    xp = []
            #    yp = []

            for iEdgeOnCell in range(0,nEdgesOnCell[iCell]):

                iEdgeOnCell1 = iEdgeOnCell
                iEdgeOnCell2 = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])

                iVertex1 = verticesOnCell[iCell,iEdgeOnCell1]
                iVertex2 = verticesOnCell[iCell,iEdgeOnCell2]

                areaTriangle = triangle_area(xCell[iCell], yCell[iCell], \
                                             xVertex[iVertex1], yVertex[iVertex1], \
                                             xVertex[iVertex2], yVertex[iVertex2])

                # subtriangle
                T11, T12, T21, T22, Cx, Cy = \
                    get_real_coords_from_barycentric_coords_coefficients(xCell[iCell], yCell[iCell], \
                                                                         xVertex[iVertex1], yVertex[iVertex1], \
                                                                         xVertex[iVertex2], yVertex[iVertex2])

                x, y = get_real_coords_from_barycentric_coords_array(u, v, T11, T12, T21, T22, Cx, Cy)
                basisWeights = np.zeros((nEdgesOnCell[iCell],nIntegrationPoints))
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    basisWeights[iVertexOnCell,:] = wachspress_basis_function_array(nEdgesOnCell[iCell],
                                                                                    iVertexOnCell,
                                                                                    x[:]-xCell[iCell],
                                                                                    y[:]-yCell[iCell],
                                                                                    wachspressKappa[:,:,iCell],
                                                                                    wachspressA[:,iCell],
                                                                                    wachspressB[:,iCell],
                                                                                    nEdgesOnCellSubset,
                                                                                    vertexIndexSubset)

                for iWeight in range(0, nIntegrationPoints):

                    e11Analytical, e22Analytical, e12Analytical, x, y = \
                        get_analytical_strain_triangle(u[iWeight], v[iWeight], T11, T12, T21, T22, Cx, Cy, xMin, yMin)

                    #if (testIntegrationPoints):
                    #    xp.append(x)
                    #    yp.append(y)

                    e11 = 0.0
                    e22 = 0.0
                    #e12 = 0.0
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                        #basisWeight = wachspress_basis_function(nEdgesOnCell[iCell],
                        #                                        iVertexOnCell,
                        #                                        x-xCell[iCell],
                        #                                        y-yCell[iCell],
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
            #    axis = fig.add_subplot(111)

            #    axis.scatter(xCell[iCell], yCell[iCell], c="red")

            #    xv = []
            #    yv = []
            #    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            #        iVertex = verticesOnCell[iCell,iVertexOnCell]
            #        axis.plot([xCell[iCell],xVertex[iVertex]],
            #                  [yCell[iCell],yVertex[iVertex]], c="g")
            #        xv.append(xVertex[iVertex])
            #        yv.append(yVertex[iVertex])
            #    iVertex = verticesOnCell[iCell,0]
            #    xv.append(xVertex[iVertex])
            #    yv.append(yVertex[iVertex])
            #    axis.plot(xv, yv)

            #    axis.scatter(xp, yp, c="k")

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

def get_norm_wachspress_integral(filenameGrid, filenameIn, strain11Name, strain22Name, strain12Name):

    fileGrid = Dataset(filenameGrid,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    areaCell = fileGrid.variables["areaCell"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]

    fileGrid.close()

    verticesOnCell[:] -= 1

    fileIn = Dataset(filenameIn,"r")

    strain11 = fileGrid.variables[strain11Name][0,:,:]
    strain22 = fileGrid.variables[strain22Name][0,:,:]
    strain12 = fileGrid.variables[strain12Name][0,:,:]

    fileIn.close()

    useCell = get_use_cell(filenameGrid)

    normE11, normE22 = L2_norm_wachspress_integral(nCells,
                                                   maxEdges,
                                                   nEdgesOnCell,
                                                   verticesOnCell,
                                                   areaCell,
                                                   xCell,
                                                   yCell,
                                                   xVertex,
                                                   yVertex,
                                                   useCell,
                                                   strain11,
                                                   strain22,
                                                   strain12)

    return normE11, normE22

#--------------------------------------------------------

def L2_norm_pwl_integral(nCells,
                         maxEdges,
                         nEdgesOnCell,
                         verticesOnCell,
                         edgesOnCell,
                         dvEdge,
                         areaCell,
                         xCell,
                         yCell,
                         xVertex,
                         yVertex,
                         useCell,
                         strain11,
                         strain22,
                         strain12):

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    xLocal, yLocal = local_coords(nCells,
                                  maxEdges,
                                  nEdgesOnCell,
                                  verticesOnCell,
                                  xCell,
                                  yCell,
                                  xVertex,
                                  yVertex)

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

    # testing
    #testIntegrationPoints = True
    #testAreaIntegration = False
    #if (testAreaIntegration):
    #    areaTest = 0.0
    #    strain11[:] = 1.0
    #    strain22[:] = 1.0
    #    strain12[:] = 1.0

    for iCell in tqdm(range(0,nCells)):

        if (useCell[iCell]):

            #if (testAreaIntegration):
            #    areaTest += areaCell[iCell]

            nEdgesOnCellSubset, vertexIndexSubset = wachspress_indexes(nEdgesOnCell[iCell])

            xC = 0.0
            yC = 0.0
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                xC += xVertex[iVertex]
                yC += yVertex[iVertex]
            xC /= float(nEdgesOnCell[iCell])
            yC /= float(nEdgesOnCell[iCell])

            #if (testIntegrationPoints):
            #    xp = []
            #    yp = []

            for iEdgeOnCell in range(0,nEdgesOnCell[iCell]):

                iEdgeOnCell1 = iEdgeOnCell
                iEdgeOnCell2 = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])

                iVertex1 = verticesOnCell[iCell,iEdgeOnCell1]
                iVertex2 = verticesOnCell[iCell,iEdgeOnCell2]

                areaTriangle = triangle_area(xC, yC, \
                                             xVertex[iVertex1], yVertex[iVertex1], \
                                             xVertex[iVertex2], yVertex[iVertex2])

                # subtriangle
                T11, T12, T21, T22, Cx, Cy = \
                    get_real_coords_from_barycentric_coords_coefficients(xC, yC, \
                                                                         xVertex[iVertex1], yVertex[iVertex1], \
                                                                         xVertex[iVertex2], yVertex[iVertex2])

                x, y = get_real_coords_from_barycentric_coords_array(u, v, T11, T12, T21, T22, Cx, Cy)
                basisWeights = np.zeros((nEdgesOnCell[iCell],nIntegrationPoints))
                for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                    iSubCell = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])
                    basisWeights[iVertexOnCell,:] = pwl_basis_function_array(nEdgesOnCell[iCell],
                                                                             iVertexOnCell,
                                                                             iSubCell,
                                                                             subBasisGradientU[iCell,:,:],
                                                                             subBasisGradientV[iCell,:,:],
                                                                             subBasisConstant[iCell,:,:],
                                                                             x[:]-xCell[iCell],
                                                                             y[:]-yCell[iCell])

                for iWeight in range(0, nIntegrationPoints):

                    e11Analytical, e22Analytical, e12Analytical, x, y = \
                        get_analytical_strain_triangle(u[iWeight], v[iWeight], T11, T12, T21, T22, Cx, Cy, xMin, yMin)

                    #if (testIntegrationPoints):
                    #    xp.append(x)
                    #    yp.append(y)

                    e11 = 0.0
                    e22 = 0.0
                    #e12 = 0.0
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                        iSubCell = wrapped_index(iEdgeOnCell+1,nEdgesOnCell[iCell])
                        #basisWeight = pwl_basis_function(nEdgesOnCell[iCell],
                        #                                 iVertexOnCell,
                        #                                 iSubCell,
                        #                                 subBasisGradientU[iCell,:,:],
                        #                                 subBasisGradientV[iCell,:,:],
                        #                                 subBasisConstant[iCell,:,:],
                        #                                 x-xCell[iCell],
                        #                                 y-yCell[iCell])
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
            #    axis = fig.add_subplot(111)

            #    axis.scatter(xCell[iCell], yCell[iCell], c="red")
            #    axis.scatter(xC, yC, c="b")

            #    xv = []
            #    yv = []
            #    for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            #        iVertex = verticesOnCell[iCell,iVertexOnCell]
            #        axis.plot([xCell[iCell],xVertex[iVertex]],
            #                  [yCell[iCell],yVertex[iVertex]], c="g")
            #        xv.append(xVertex[iVertex])
            #        yv.append(yVertex[iVertex])
            #    iVertex = verticesOnCell[iCell,0]
            #    xv.append(xVertex[iVertex])
            #    yv.append(yVertex[iVertex])
            #    axis.plot(xv, yv)

            #    axis.scatter(xp, yp, c="k")

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

def get_norm_pwl_integral(filenameGrid, filenameIn, strain11Name, strain22Name, strain12Name):

    fileGrid = Dataset(filenameGrid,"r")

    nCells = len(fileGrid.dimensions["nCells"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])
    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    edgesOnCell = fileGrid.variables["edgesOnCell"][:]
    dvEdge = fileGrid.variables["dvEdge"][:]
    areaCell = fileGrid.variables["areaCell"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]

    fileGrid.close()

    verticesOnCell[:] -= 1
    edgesOnCell[:] -= 1

    fileIn = Dataset(filenameIn,"r")

    strain11 = fileGrid.variables[strain11Name][0,:,:]
    strain22 = fileGrid.variables[strain22Name][0,:,:]
    strain12 = fileGrid.variables[strain12Name][0,:,:]

    fileIn.close()

    useCell = get_use_cell(filenameGrid)

    normE11, normE22 = L2_norm_pwl_integral(nCells,
                                            maxEdges,
                                            nEdgesOnCell,
                                            verticesOnCell,
                                            edgesOnCell,
                                            dvEdge,
                                            areaCell,
                                            xCell,
                                            yCell,
                                            xVertex,
                                            yVertex,
                                            useCell,
                                            strain11,
                                            strain22,
                                            strain12)

    return normE11, normE22

#--------------------------------------------------------

def L2_norm_weak_integral(nVertices,
                          vertexDegree,
                          edgesOnVertex,
                          cellsOnEdge,
                          cellsOnVertex,
                          dcEdge,
                          useVertex,
                          xCell,
                          yCell,
                          xVertex,
                          yVertex,
                          strain11,
                          strain22,
                          strain12):

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    integrationOrder = 7
    nIntegrationPoints, u, weights = integration_weights_gauss_lobatto(integrationOrder)

    normE11  = 0.0
    denomE11 = 0.0
    normE22  = 0.0
    denomE22 = 0.0
    #normE12  = 0.0
    #denomE12 = 0.0

    #testIntegrationPoints = True
    #testAreaIntegration = True
    #if (testAreaIntegration):
    #    lengthTest = 0.0
    #    strain11[:] = 1.0
    #    strain22[:] = 1.0
    #    strain12[:] = 1.0

    for iVertex in tqdm(range(0,nVertices)):

        if (useVertex[iVertex]):

            #if (testIntegrationPoints):
            #    xp = []
            #    yp = []

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

                    x = w1 * xCell[iCell1] + w2 * xCell[iCell2]
                    y = w1 * yCell[iCell1] + w2 * yCell[iCell2]

                    #if (testIntegrationPoints):
                    #    xp.append(x)
                    #    yp.append(y)

                    uVel, vVel, e11Analytical, e22Analytical, e12Analytical, divu, divv = velocities_strains_stress_divergences(x-xMin, y-yMin, "sinusoid")

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
            #    axis = fig.add_subplot(111)

            #    axis.scatter(xVertex[iVertex], yVertex[iVertex], c="red")

            #    xv = []
            #    yv = []
            #    for iCellOnVertex in range(0,vertexDegree):
            #        iCell = cellsOnVertex[iVertex,iCellOnVertex]
            #        axis.plot([xVertex[iVertex],xCell[iCell]],
            #                  [yVertex[iVertex],yCell[iCell]], c="g")
            #        xv.append(xCell[iCell])
            #        yv.append(yCell[iCell])
            #    iCell = cellsOnVertex[iVertex,0]
            #    xv.append(xCell[iCell])
            #    yv.append(yCell[iCell])
            #    axis.plot(xv, yv, c="b")

            #    axis.scatter(xp, yp, c="k")

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

def get_norm_weak_integral(filenameGrid, filenameIn, strain11Name, strain22Name, strain12Name):

    fileGrid = Dataset(filenameGrid,"r")

    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    edgesOnVertex = fileGrid.variables["edgesOnVertex"][:]
    cellsOnEdge = fileGrid.variables["cellsOnEdge"][:]
    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    dcEdge = fileGrid.variables["dcEdge"][:]
    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]

    fileGrid.close()

    edgesOnVertex[:] -= 1
    cellsOnEdge[:] -= 1
    cellsOnVertex[:] -= 1

    fileIn = Dataset(filenameIn,"r")

    strain11 = fileGrid.variables[strain11Name][0,:]
    strain22 = fileGrid.variables[strain22Name][0,:]
    strain12 = fileGrid.variables[strain12Name][0,:]

    fileIn.close()

    useVertex = get_use_vertex(filenameGrid)

    normE11, normE22 = L2_norm_weak_integral(nVertices,
                                             vertexDegree,
                                             edgesOnVertex,
                                             cellsOnEdge,
                                             cellsOnVertex,
                                             dcEdge,
                                             useVertex,
                                             xCell,
                                             yCell,
                                             xVertex,
                                             yVertex,
                                             strain11,
                                             strain22,
                                             strain12)

    return normE11, normE22

#--------------------------------------------------------

def L2_norm_weak(numerical, analytical, nCells, areaCell, useCell):

    degreesToRadians = pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (useCell[iCell]):

            norm  = norm  + areaCell[iCell] * pow(numerical[iCell] - analytical[iCell],2)

            denom = denom + areaCell[iCell] * pow(analytical[iCell],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_weak(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    areaCell = fileMPAS.variables["areaCell"][:]

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    strain22weak = fileMPAS.variables["strain22weak"][0,:]
    strain12weak = fileMPAS.variables["strain12weak"][0,:]

    useCell = get_use_cell(filename)

    normdfdx = L2_norm_weak(strain11weak, strain11CellAnalytical, nCells, areaCell, useCell)
    normdfdy = L2_norm_weak(strain22weak, strain22CellAnalytical, nCells, areaCell, useCell)

    fileMPAS.close()

    return normdfdx, normdfdy

#--------------------------------------------------------

def L2_norm_var(numerical, analytical, nCells, nEdgesOnCell, verticesOnCell, areaCell, useCell):

    degreesToRadians = pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iCell in range(0,nCells):

        if (useCell[iCell]):

            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

                iVertex = verticesOnCell[iCell,iVertexOnCell]

                norm  = norm  + areaCell[iCell] * pow(numerical[iCell,iVertexOnCell] - analytical[iVertex],2)

                denom = denom + areaCell[iCell] * pow(analytical[iVertex],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_var(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    areaCell = fileMPAS.variables["areaCell"][:]

    nEdgesOnCell = fileMPAS.variables["nEdgesOnCell"][:]
    verticesOnCell = fileMPAS.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    strain11var = fileMPAS.variables["strain11var"][0,:]
    strain22var = fileMPAS.variables["strain22var"][0,:]
    strain12var = fileMPAS.variables["strain12var"][0,:]

    useCell = get_use_cell(filename)

    normdfdx = L2_norm_var(strain11var, strain11VertexAnalytical, nCells, nEdgesOnCell, verticesOnCell, areaCell, useCell)
    normdfdy = L2_norm_var(strain22var, strain22VertexAnalytical, nCells, nEdgesOnCell, verticesOnCell, areaCell, useCell)

    fileMPAS.close()

    return normdfdx, normdfdy

#--------------------------------------------------------

def L2_norm_var_avg(numerical, analytical, nVertices, areaTriangle, useVertex):

    degreesToRadians = pi / 180.0

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (useVertex[iVertex]):

            norm  = norm  + areaTriangle[iVertex] * pow(numerical[iVertex] - analytical[iVertex],2)

            denom = denom + areaTriangle[iVertex] * pow(analytical[iVertex],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_var_avg(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])

    areaTriangle = fileMPAS.variables["areaTriangle"][:]

    strain11varAvg = fileMPAS.variables["strain11varAvg"][0,:]
    strain22varAvg = fileMPAS.variables["strain22varAvg"][0,:]
    strain12varAvg = fileMPAS.variables["strain12varAvg"][0,:]

    useVertex = get_use_vertex(filename)

    normdfdx = L2_norm_var_avg(strain11varAvg, strain11VertexAnalytical, nVertices, areaTriangle, useVertex)
    normdfdy = L2_norm_var_avg(strain22varAvg, strain22VertexAnalytical, nVertices, areaTriangle, useVertex)

    fileMPAS.close()

    return normdfdx, normdfdy

#--------------------------------------------------------

def L2_norm_weak_avg(numerical, analytical, nVertices, areaTriangle, useVertex):

    norm  = 0.0
    denom = 0.0

    for iVertex in range(0,nVertices):

        if (useVertex[iVertex]):

            norm  += areaTriangle[iVertex] * pow(numerical[iVertex] - analytical[iVertex],2)

            denom += areaTriangle[iVertex] * pow(analytical[iVertex],2)

    norm = sqrt(norm / denom)

    return norm

#--------------------------------------------------------

def get_norm_weak_avg(filenameIC, filename):

    fileIC = Dataset(filenameIC, "r")

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    fileIC.close()

    fileMPAS = Dataset(filename, "r")

    nVertices = len(fileMPAS.dimensions["nVertices"])
    nCells = len(fileMPAS.dimensions["nCells"])
    vertexDegree = len(fileMPAS.dimensions["vertexDegree"])

    areaTriangle = fileMPAS.variables["areaTriangle"][:]
    areaCell = fileMPAS.variables["areaCell"][:]

    cellsOnVertex = fileMPAS.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    strain11weak = fileMPAS.variables["strain11weak"][0,:]
    strain22weak = fileMPAS.variables["strain22weak"][0,:]
    strain12weak = fileMPAS.variables["strain12weak"][0,:]

    strain11weakAvg = np.zeros(nVertices)
    strain22weakAvg = np.zeros(nVertices)
    strain12weakAvg = np.zeros(nVertices)

    for iVertex in range(0,nVertices):

        strain11 = 0.0
        strain22 = 0.0
        strain12 = 0.0
        denominator = 0.0

        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell < nCells):
                strain11 += areaCell[iCell] * strain11weak[iCell]
                strain22 += areaCell[iCell] * strain22weak[iCell]
                strain12 += areaCell[iCell] * strain12weak[iCell]
                denominator += areaCell[iCell]

        if (denominator > 0.0):
            strain11weakAvg[iVertex] = strain11 / denominator
            strain22weakAvg[iVertex] = strain22 / denominator
            strain12weakAvg[iVertex] = strain12 / denominator

    useVertex = get_use_vertex(filename)

    normdfdx = L2_norm_var_avg(strain11weakAvg, strain11VertexAnalytical, nVertices, areaTriangle, useVertex)
    normdfdy = L2_norm_var_avg(strain22weakAvg, strain22VertexAnalytical, nVertices, areaTriangle, useVertex)

    fileMPAS.close()

    return normdfdx, normdfdy

#--------------------------------------------------------

def get_resolution(filename):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    nEdges = len(fileMPAS.dimensions["nEdges"])

    degreesToRadians = pi / 180.0

    cellsOnEdge = fileMPAS.variables["cellsOnEdge"][:]
    cellsOnEdge[:] = cellsOnEdge[:] - 1

    dcEdge = fileMPAS.variables["dcEdge"][:]
    latEdge = fileMPAS.variables["latEdge"][:]

    resolution = 0.0
    denom = 0.0
    for iEdge in range(0,nEdges):
        iCell1 = cellsOnEdge[iEdge,0]
        iCell2 = cellsOnEdge[iEdge,1]
        if (iCell1 != -1 and iCell2 != -1):
            resolution = resolution + dcEdge[iEdge]
            denom = denom + 1.0

    resolution = resolution / denom

    fileMPAS.close()

    return resolution

#--------------------------------------------------------

def get_use_cell(filenameIn):

    fileIn = Dataset(filenameIn,"r")
    nCells = len(fileIn.dimensions["nCells"])
    nVertices = len(fileIn.dimensions["nVertices"])
    nEdgesOnCell = fileIn.variables["nEdgesOnCell"][:]
    cellsOnCell = fileIn.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1
    verticesOnCell = fileIn.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1
    interiorCell = fileIn.variables["interiorCell"][0,:]
    fileIn.close()

    useCell = np.ones(nCells,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            useCell[iCell] = 0
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                useCell[iCell2] = 0

    return useCell

#--------------------------------------------------------

def get_use_vertex(filenameIn):

    fileIn = Dataset(filenameIn,"r")
    nCells = len(fileIn.dimensions["nCells"])
    nVertices = len(fileIn.dimensions["nVertices"])
    nEdgesOnCell = fileIn.variables["nEdgesOnCell"][:]
    cellsOnCell = fileIn.variables["cellsOnCell"][:]
    cellsOnCell[:] = cellsOnCell[:] - 1
    verticesOnCell = fileIn.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1
    interiorCell = fileIn.variables["interiorCell"][0,:]
    fileIn.close()

    useVertex = np.ones(nVertices,dtype="i")
    for iCell in range(0,nCells):
        if (interiorCell[iCell] == 0):
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                useVertex[iVertex] = 0
            for iCellOnCell in range(0,nEdgesOnCell[iCell]):
                iCell2 = cellsOnCell[iCell,iCellOnCell]
                if (iCell2 < nCells):
                    for iVertexOnCell in range(0,nEdgesOnCell[iCell2]):
                        iVertex = verticesOnCell[iCell2,iVertexOnCell]
                        useVertex[iVertex] = 0

    return useVertex

#--------------------------------------------------------

def strain_scaling():

    strains = ["dfdx","dfdy"]
    iStrainDfdx = 0
    iStrainDfdy = 1

    gridTypes = ["quad","hex"]
    #gridTypes = ["hex"]



    operatorMethods = ["wachspress","pwl","weak","wachspress_avg","pwl_avg","weak_avg"]
    #operatorMethods = ["weak"]

    grids = {"hex" :["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}
    #grids = {"hex" :["0082x0094"],
    #         "quad":["0080x0080"]}


    dirname = {"wachspress":"wachspress",
               "pwl":"pwl",
               "weak":"weak",
               "wachspress_avg":"wachspress",
               "pwl_avg":"pwl",
               "weak_avg":"weak"}

    legendLabels = {"wachspress":"Wachspress",
                    "pwl":"PWL",
                    "weak":"FV",
                    "wachspress_avg":"Wachs. Avg",
                    "pwl_avg":"PWL Avg",
                    "weak_avg":"FV Avg"}

    lineColours = {"wachspress":"black",
                   "pwl":"grey",
                   "weak":"red",
                   "wachspress_avg":"blue",
                   "pwl_avg":"green",
                   "weak_avg":"darkturquoise"}

    markers = {"wachspress":"+",
               "pwl":"x",
               "weak":"^",
               "wachspress_avg":"+",
               "pwl_avg":"x",
               "weak_avg":"^"}

    #lineStyles = {"hex":"solid",
    #              "quad":"dashed"}
    lineStyles = {"hex":"solid",
                  "quad":"solid"}

    strainLabels = {"dfdx":{"hex":r"(c) Hex. $\partial{f}/\partial{x}$",
                            "quad":r"(a) Quad. $\partial{f}/\partial{x}$"},
                    "dfdy":{"hex":r"(d) Hex. $\partial{f}/\partial{y}$",
                            "quad":r"(b) Quad. $\partial{f}/\partial{y}$"}}

    xlabels = {"hex":"Grid resolution",
               "quad":None}

    ylabels = {"dfdx":r"$L_2$ error norm",
               "dfdy":None}


    cm = 1/2.54  # centimeters in inches
    plt.rc('font', family="Times New Roman")
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


    fig, axes = plt.subplots(2,2,figsize=(15*cm,14*cm))

    #xMin = 6e-3
    #xMax = 1e-2
    #yMin = 4e-4
    xMin = 1.6e-3
    xMax = 3.0e-3
    yMin = 8.0e-4

    # linear scaling
    scale = yMin / pow(xMin,1)
    scaleMinLin = pow(xMin,1) * scale
    scaleMaxLin = pow(xMax,1) * scale

    # quadratic scaling
    scale = yMin / pow(xMin,2)
    scaleMinQuad = pow(xMin,2) * scale
    scaleMaxQuad = pow(xMax,2) * scale

    plotHandles = []

    iGrid = 0
    for gridType in gridTypes:

        print("   Grid type: ", iGrid, gridType)

        if ("dfdx" in strains):
            axes[iGrid,iStrainDfdx].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
            axes[iGrid,iStrainDfdx].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')
        if ("dfdy" in strains):
            axes[iGrid,iStrainDfdy].loglog([xMin, xMax], [scaleMinLin,scaleMaxLin], linestyle=':', color='k')
            axes[iGrid,iStrainDfdy].loglog([xMin, xMax], [scaleMinQuad,scaleMaxQuad], linestyle=':', color='k')

        iPlot = 0
        for operatorMethod in operatorMethods:

            print("     Operator Method: ", operatorMethod)

            x = []
            ydfdx = []
            ydfdy = []

            for grid in grids[gridType]:

                print("       Grid: ", grid)

                filenameIC = "ic_%s_%s.nc" %(gridType,grid)
                filename = "output_%s_%s_%s/output.2000.nc" %(gridType, dirname[operatorMethod], grid)

                print("         ", filename, filenameIC)

                if (operatorMethod == "wachspress"):
                    #normdfdx, normdfdy = get_norm_var(filenameIC, filename)
                    normdfdx, normdfdy = get_norm_wachspress_integral(filename, filename, "strain11var", "strain22var", "strain12var")
                elif (operatorMethod == "pwl"):
                    #normdfdx, normdfdy = get_norm_var(filenameIC, filename)
                    normdfdx, normdfdy = get_norm_pwl_integral(filename, filename, "strain11var", "strain22var", "strain12var")
                elif (operatorMethod == "weak"):
                    #normdfdx, normdfdy = get_norm_weak(filenameIC, filename)
                    normdfdx, normdfdy = get_norm_weak_integral(filename, filename, "strain11weak", "strain22weak", "strain12weak")
                elif (operatorMethod == "wachspress_avg"):
                    #normdfdx, normdfdy = get_norm_var_avg(filenameIC, filename)
                    normdfdx, normdfdy = get_norm_wachspress_integral(filename, filename, "strain11varAvg", "strain22varAvg", "strain12varAvg")
                elif (operatorMethod == "pwl_avg"):
                    #normdfdx, normdfdy = get_norm_var_avg(filenameIC, filename)
                    normdfdx, normdfdy = get_norm_pwl_integral(filename, filename, "strain11varAvg", "strain22varAvg", "strain12varAvg")
                elif (operatorMethod == "weak_avg"):
                    #normdfdx, normdfdy = get_norm_weak_avg(filenameIC, filename)
                    normdfdx, normdfdy = get_norm_wachspress_integral(filename, filename, "strain11weakAvg", "strain22weakAvg", "strain12weakAvg")

                x.append(get_resolution(filename))
                if ("dfdx" in strains):
                    ydfdx.append(normdfdx)
                if ("dfdy" in strains):
                    ydfdy.append(normdfdy)

            #if (gridType == "hex"):
            legendLabel = "%s" %(legendLabels[operatorMethod])
            #else:
            #    legendLabel = "_nolegend_"
            if ("dfdx" in strains):
                plotHandle, = axes[iGrid,iStrainDfdx].loglog(x, ydfdx, marker=markers[operatorMethod], color=lineColours[operatorMethod], ls=lineStyles[gridType], markersize=5.0, label=legendLabel, fillstyle=None)
            if ("dfdy" in strains):
                plotHandle, = axes[iGrid,iStrainDfdy].loglog(x, ydfdy, marker=markers[operatorMethod], color=lineColours[operatorMethod], ls=lineStyles[gridType], markersize=5.0, label=legendLabel, fillstyle=None)
            plotHandles.append(plotHandle)

            iPlot = iPlot + 1

        # dfdx
        if ("dfdx" in strains):
            legend1 = axes[iGrid,iStrainDfdx].legend(handles=plotHandles[0:3], frameon=False, loc=2, fontsize=8, handlelength=2)
            legend2 = axes[iGrid,iStrainDfdx].legend(handles=plotHandles[3:6], frameon=False, loc=4, fontsize=8, handlelength=2)
            axes[iGrid,iStrainDfdx].add_artist(legend1)
            #axes[iGrid,iStrainDfdx].legend(frameon=False, loc=4, fontsize=8, handlelength=2)

            axes[iGrid,iStrainDfdx].set_xlabel(xlabels[gridType])
            axes[iGrid,iStrainDfdx].set_ylabel(ylabels[strains[iStrainDfdx]])
            axes[iGrid,iStrainDfdx].set_title(strainLabels[strains[iStrainDfdx]][gridType],loc="left")
            axes[iGrid,iStrainDfdx].set_ylim(None,None)
            axes[iGrid,iStrainDfdx].set_xticks(ticks=[3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3],minor=True)
            axes[iGrid,iStrainDfdx].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
            axes[iGrid,iStrainDfdx].set_xticks(ticks=[2e-3,1e-2],minor=False)
            axes[iGrid,iStrainDfdx].set_xticklabels(labels=[r'$2\times 10^{-3}$',r'$10^{-2}$'],minor=False)

        # dfdy
        if ("dfdy" in strains):
            legend1 = axes[iGrid,iStrainDfdy].legend(handles=plotHandles[0:3], frameon=False, loc=2, fontsize=8, handlelength=2)
            legend2 = axes[iGrid,iStrainDfdy].legend(handles=plotHandles[3:6], frameon=False, loc=4, fontsize=8, handlelength=2)
            axes[iGrid,iStrainDfdy].add_artist(legend1)
            #axes[iGrid,iStrainDfdy].legend(frameon=False, loc=4, fontsize=8, handlelength=2)

            axes[iGrid,iStrainDfdy].set_xlabel(xlabels[gridType])
            axes[iGrid,iStrainDfdy].set_ylabel(ylabels[strains[iStrainDfdy]])
            axes[iGrid,iStrainDfdy].set_title(strainLabels[strains[iStrainDfdy]][gridType],loc="left")
            axes[iGrid,iStrainDfdy].set_ylim(None,None)
            axes[iGrid,iStrainDfdy].set_xticks(ticks=[3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3],minor=True)
            axes[iGrid,iStrainDfdy].set_xticklabels(labels=[None,None,None,None,None,None,None],minor=True)
            axes[iGrid,iStrainDfdy].set_xticks(ticks=[2e-3,1e-2],minor=False)
            axes[iGrid,iStrainDfdy].set_xticklabels(labels=[r'$2\times 10^{-3}$',r'$10^{-2}$'],minor=False)

        iGrid += 1

    plt.tight_layout(pad=0.2, w_pad=0.6, h_pad=0.2)
    plt.savefig("derivative_scaling.png", dpi=300)
    plt.savefig("derivative_scaling.eps")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_scaling()
