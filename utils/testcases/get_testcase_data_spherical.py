import os

#-------------------------------------------------------------------------------

def get_testcase_data_spherical(getGridFiles = True,
                                getGraphFiles = True):

    MPAS_SEAICE_STANDALONE_DATA = os.environ.get('MPAS_SEAICE_STANDALONE_DATA')
    if (MPAS_SEAICE_STANDALONE_DATA is None):
        raise Exception("MPAS_SEAICE_STANDALONE_DATA must be set")

    nCellsArray = [10242, 163842, 2562, 40962]

    for nCells in nCellsArray:

        # grid file
        if (getGridFiles):
            filanameDst = "grid.%i.nc" %(nCells)
            filenameSrc = "%s/testcases/strain_stress_divergence/%s" %(MPAS_SEAICE_STANDALONE_DATA,filanameDst)

            if (not os.path.isfile(filanameDst)):
                os.symlink(filenameSrc, filanameDst)

        # graph file
        if (getGraphFiles):
            filanameDst = "graph.%i.info" %(nCells)
            filenameSrc = "%s/testcases/strain_stress_divergence/%s" %(MPAS_SEAICE_STANDALONE_DATA,filanameDst)

            if (not os.path.isfile(filanameDst)):
                os.symlink(filenameSrc, filanameDst)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data_spherical()
