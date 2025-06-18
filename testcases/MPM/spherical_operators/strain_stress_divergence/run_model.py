import sys

sys.path.append("../../../utils/testcases")
from log_messages import log_message
from execute_model import execute_model

import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model(logFile=None):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    gridSizes = [2562, 10242, 40962, 163842]
    operatorMethods = ["wachspress","pwl","weak","weakwachs","weakpwl","wachspress_alt","pwl_alt","weakwachs_alt","mpmvar","mpmweak"]
    #operatorMethods = ["mpmvar","mpmweak"]

    for operatorMethod in operatorMethods:

        print("Operator Method: ", operatorMethod)

        for gridSize in gridSizes:

            print("  Gridsize: ", gridSize)

            os.system("rm grid.nc ic.nc particles.nc")
            os.system("ln -s grid.%i.nc grid.nc" %(gridSize))
            os.system("ln -s ic_%i.nc ic.nc" %(gridSize))
            os.system("ln -s particles_%i.nc particles.nc" %(gridSize))

            if (operatorMethod == "wachspress"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress",
                                                "config_variational_denominator_type":"original"}}
            elif (operatorMethod == "pwl"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"pwl",
                                                "config_variational_denominator_type":"original"}}
            elif (operatorMethod == "weak"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                "config_stress_divergence_scheme":"weak",
                                                "config_variational_denominator_type":"original"}}
            elif (operatorMethod == "weakwachs"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress",
                                                "config_variational_denominator_type":"original"}}
            elif (operatorMethod == "weakpwl"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"pwl",
                                                "config_variational_denominator_type":"original"}}
            elif (operatorMethod == "wachspress_alt"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress",
                                                "config_variational_denominator_type":"alternate"}}
            elif (operatorMethod == "pwl_alt"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"pwl",
                                                "config_variational_denominator_type":"alternate"}}
            elif (operatorMethod == "weakwachs_alt"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                "config_stress_divergence_scheme":"variational",
                                                "config_variational_basis":"wachspress",
                                                "config_variational_denominator_type":"alternate"}}
            elif (operatorMethod == "mpmvar"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                "config_stress_divergence_scheme":"variational"},
                            "use_sections": {"config_use_mpm": True}}
            elif (operatorMethod == "mpmweak"):
                nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                "config_stress_divergence_scheme":"weak"},
                            "use_sections": {"config_use_mpm": True}}

            f90nml.patch("namelist.seaice.strain_stress_divergence", nmlPatch, "namelist.seaice.%s.%i" %(operatorMethod, gridSize))

            os.system("rm -rf namelist.seaice streams.seaice output_%s_%i" %(operatorMethod, gridSize))
            os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(operatorMethod, gridSize))
            os.system("ln -s streams.seaice.strain_stress_divergence streams.seaice")

            execute_model(MPAS_SEAICE_EXECUTABLE,
                          1,
                          logFile)

            os.system("mv output output_%s_%i" %(operatorMethod, gridSize))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
