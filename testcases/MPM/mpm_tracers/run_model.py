import os
import sys
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

sys.path.append("../../../utils/testcases")
from log_messages import log_message
from execute_model import execute_model

sys.path.append("../../../testing")
from testing_utils import add_pio_namelist_changes, create_new_namelist

#-------------------------------------------------------------------------------

def run_model(logFileOverview):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    MPAS_SEAICE_METIS_PATH = os.environ.get('MPAS_SEAICE_METIS_PATH')
    if (MPAS_SEAICE_METIS_PATH is None):
        MPAS_SEAICE_METIS_PATH = "gpmetis"

    operatorMethods = ["mpm_tracers", "mpas"]

    for operatorMethod in operatorMethods:
        print("   operatorMethod: ", operatorMethod)

        if (operatorMethod == "mpm_tracers"):
            nmlPatch = {"use_sections": {"config_use_mpm": True},
                                 "mpm": {"config_use_mpm_tracers": True},
                      "column_package": {"config_column_element_type":"particles"},
                           "advection": {"config_advection_type":"mpm"}}
        elif (operatorMethod == "mpas"):
            nmlPatch = {"use_sections": {"config_use_mpm": False},
                                 "mpm": {"config_use_mpm_tracers": False},
                      "column_package": {"config_column_element_type":"cells"},
                           "advection": {"config_advection_type":"incremental_remap"}}

        f90nml.patch("namelist.seaice.default", nmlPatch, "namelist.seaice.%s" %(operatorMethod))

        os.system("rm -rf namelist.seaice")
        os.system("ln -s namelist.seaice.%s namelist.seaice" %(operatorMethod))

        if (not os.path.isdir("output")):
            os.mkdir("output")

        cmd = ("rm -rf output_%s" %(operatorMethod))
        print(cmd)
        os.system(cmd)

        execute_model(MPAS_SEAICE_EXECUTABLE,
                      1,
                      logFileOverview)

        cmd = "mv output output_%s" %(operatorMethod)
        print(cmd)
        os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
