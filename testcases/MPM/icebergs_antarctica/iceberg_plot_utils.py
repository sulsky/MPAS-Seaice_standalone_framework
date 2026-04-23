
#-------------------------------------------------------------------------------

def projection_scalar(x, y, z, location):

    if (location == "antarctica"):

        xPlot = y
        yPlot = x

    elif (location == "greenland"):

        xPlot = x - 1.25e6
        yPlot = y + 1e6

    return xPlot, yPlot

#-------------------------------------------------------------------------------

def projection_list(x, y, z, location):

    if (location == "antarctica"):

        xPlot = y
        yPlot = x

    elif (location == "greenland"):

        xPlot = [xi - 1.25e6 for xi in x]
        yPlot = [yi + 1.e6   for yi in y]

    return xPlot, yPlot

#-------------------------------------------------------------------------------

def plot_limits(xMin, xMax, yMin, yMax, enlarge=1.05):

    dx = xMax-xMin
    dy = yMax-yMin
    dxy = max(dx,dy) * 0.5
    xc = 0.5*(xMin+xMax)
    yc = 0.5*(yMin+yMax)
    xMin = xc - dxy*enlarge
    xMax = xc + dxy*enlarge
    yMin = yc - dxy*enlarge
    yMax = yc + dxy*enlarge

    return xMin, xMax, yMin, yMax

#-------------------------------------------------------------------------------
