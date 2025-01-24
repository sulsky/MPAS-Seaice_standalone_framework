import os
import sys
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

sys.path.append("../../../testing")
from testing_utils import add_pio_namelist_changes, create_new_namelist

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    MPAS_SEAICE_METIS_PATH = os.environ.get('MPAS_SEAICE_METIS_PATH')
    if (MPAS_SEAICE_METIS_PATH is None):
        MPAS_SEAICE_METIS_PATH = "gpmetis"

    operatorMethods = ["mpmvar", "mpas"]

    for operatorMethod in operatorMethods:
       print("   operatorMethod: ", operatorMethod)

       if (operatorMethod == "mpmvar"):
          nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                          "config_stress_divergence_scheme":"variational"},
                      "use_sections": {"config_use_mpm": True}}
       elif (operatorMethod == "mpas"):
          nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                          "config_stress_divergence_scheme":"variational"},
                      "use_sections": {"config_use_mpm": False}}

       f90nml.patch("namelist.seaice.default", nmlPatch, "namelist.seaice.%s" %(operatorMethod))

       os.system("rm -rf namelist.seaice")
       os.system("ln -s namelist.seaice.%s namelist.seaice" %(operatorMethod))


       if (not os.path.isdir("output")):
              os.mkdir("output")

       cmd = ("rm -rf output_%s" %(operatorMethod))
       print(cmd)
       os.system(cmd)

       cmd = "mpirun -np 1 %s" %(MPAS_SEAICE_EXECUTABLE)
       #cmd = "srun --nodes=1 --cpus-per-task=1 --ntasks-per-node=%i %s" %(nProc, MPAS_SEAICE_EXECUTABLE)
       print(cmd)
       os.system(cmd)

       cmd = "mv output output_%s" %(operatorMethod)
       print(cmd)
       os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
