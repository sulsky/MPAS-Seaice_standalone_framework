import sys

sys.path.append("../../../utils/MPM/particle_initialization/")
from empty_particle_file import empty_particle_file

sys.path.append("../../../utils/MPM/icebergs/")
from create_mpas_antarctic_calving_rates import create_mpas_antarctic_calving_rates

#-------------------------------------------------------------------------------

def run_testcase():

    # create the empty iceberg file
    empty_particle_file("icebergs.nc",
                        "nIcebergs")

    # create the calving rate input file
    create_mpas_antarctic_calving_rates("grid.nc",
                                        "calving_mpas.nc")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
