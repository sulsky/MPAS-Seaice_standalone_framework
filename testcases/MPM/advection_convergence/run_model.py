import sys

sys.path.append("../../../utils/testcases")
from execute_model import execute_model

import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model(usePolympo,
              logFile=None):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    icTypes = ["cosine_bell","slotted_cylinder"]
    #icTypes = ["cosine_bell"]

    gridSizes = [2562, 10242, 40962, 163842]
    #gridSizes = [2562]

    if (usePolympo):
         usePolympoStr = "polympo"
    else:
         usePolympoStr = "nonpolympo"

    for icType in icTypes:

         for gridSize in gridSizes:

             print()
             message = "IC type: %s, Gridsize: %i" %(icType, gridSize)
             print(message)
             print("-"*len(message))

             if (not os.path.isdir("output")):
                 os.mkdir("output")

             os.system("rm grid.nc ic.nc particles.nc namelist.seaice")
             os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(usePolympoStr, gridSize))
             os.system("ln -s grid.%i.nc grid.nc" %(gridSize))
             os.system("ln -s ic_%s_%i.nc ic.nc" %(icType, gridSize))
             os.system("ln -s particles_%s_%i.nc particles.nc" %(icType, gridSize))

             os.system("rm -rf output_%s_%i_%s" %(icType, gridSize, usePolympoStr))

             execute_model(MPAS_SEAICE_EXECUTABLE,
                           1,
                           logFile)

             os.system("mv output output_%s_%i_%s" %(icType, gridSize, usePolympoStr))
             os.system("mv log.seaice.0000.out log.seaice.0000.out_%s_%i_%s" %(icType, gridSize, usePolympoStr))

#-----------------------------------------------------------------------------

if __name__ == "__main__":

    run_model(False)
