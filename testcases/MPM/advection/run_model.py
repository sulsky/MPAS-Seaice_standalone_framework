import os
import argparse

#-------------------------------------------------------------------------------

def run_model(nCells):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    if (not os.path.isdir("output")):
        os.mkdir("output")

    cmd = "rm grid.nc ic.nc"
    print(cmd)
    os.system(cmd)

    cmd = "ln -s grid.%i.nc grid.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "ln -s ic_slotted_cylinder_%i.nc ic.nc" %(nCells)
    print(cmd)
    os.system(cmd)

    cmd = "%s" %(MPAS_SEAICE_EXECUTABLE)
    print(cmd)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest="nCells", type=int)
    args = parser.parse_args()

    run_model(args.nCells)
