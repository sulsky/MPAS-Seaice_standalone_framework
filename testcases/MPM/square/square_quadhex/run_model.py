import os
import sys

try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

sys.path.append("../../../../utils/testcases")
from execute_model import execute_model

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        MPAS_SEAICE_EXECUTABLE = "../../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        print("Using executable in standard location: %s" %(MPAS_SEAICE_EXECUTABLE))

    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    operatorMethods = ["wachspress","weak","mpmvar","mpmweak","mpmmpm"]
    #operatorMethods = ["wachspress","pwl","weak","wachsavg","pwlavg","weakwachs","weakpwl"]
    #operatorMethods = ["none"]

    gridTypes = ["hex","quad"]
    #gridTypes = ["quad"]

    #grids = {"hex" :["0082x0094",
    #                 "0164x0188",
    #                 "0328x0376",
    #                 "0656x0752"],
    #         "quad":["0080x0080",
    #                 "0160x0160",
    #                 "0320x0320",
    #                 "0640x0640"]}
    grids = {"hex" :["0082x0094"],
             "quad":["0080x0080"]}
    #grids = {"quad":["0080x0080"]}


    #subcycleNumbers = [120,240,480,960,1920,3840,7680,15360,30720]
    #subcycleNumbers = [120,240,480,960,1920,3840,7680]
    subcycleNumbers = [120]

    for gridType in gridTypes:

        print("Grid type: ", gridType)

        for operatorMethod in operatorMethods:

            print("  Operator Method: ", operatorMethod)

            for grid in grids[gridType]:

                print("    Grid: ", grid)

                os.system("rm grid.nc")
                os.system("rm ic.nc")
                os.system("rm particles.nc")
                os.system("ln -s grid_%s_%s.nc grid.nc" %(gridType,grid))
                os.system("ln -s ic_%s_%s.nc ic.nc" %(gridType,grid))
                os.system("ln -s particles_%s_%s.nc particles.nc" %(gridType,grid))

                for subcycleNumber in subcycleNumbers:

                    print("      Subcycle number: ", subcycleNumber)

                    if (operatorMethod == "wachspress"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber}}
                    elif (operatorMethod == "pwl"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"pwl",
                                                        "config_elastic_subcycle_number":subcycleNumber}}
                    elif (operatorMethod == "weak"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                        "config_stress_divergence_scheme":"weak",
                                                        "config_elastic_subcycle_number":subcycleNumber}}
                    elif (operatorMethod == "wachsavg"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber,
                                                        "config_average_variational_strain":True}}
                    elif (operatorMethod == "pwlavg"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"pwl",
                                                        "config_elastic_subcycle_number":subcycleNumber,
                                                        "config_average_variational_strain":True}}
                    elif (operatorMethod == "weakwachs"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber}}
                    elif (operatorMethod == "weakpwl"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"pwl",
                                                        "config_elastic_subcycle_number":subcycleNumber}}
                    elif (operatorMethod == "mpmvar"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber},
                                    "use_sections":{"config_use_mpm":True},
                                    "column_package":{"config_column_element_type":"particles"}}
                    elif (operatorMethod == "mpmweak"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                        "config_stress_divergence_scheme":"weak",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber},
                                    "use_sections":{"config_use_mpm":True},
                                    "column_package":{"config_column_element_type":"particles"}}
                    elif (operatorMethod == "mpmmpm"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                        "config_stress_divergence_scheme":"mpm",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber},
                                    "use_sections":{"config_use_mpm":True},
                                    "column_package":{"config_column_element_type":"particles"}}
                    elif (operatorMethod == "none"):
                        nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                        "config_stress_divergence_scheme":"variational",
                                                        "config_variational_basis":"wachspress",
                                                        "config_elastic_subcycle_number":subcycleNumber}}


                    f90nml.patch("namelist.seaice.square", nmlPatch, "namelist.seaice.%s.%i" %(operatorMethod, subcycleNumber))

                    os.system("rm -rf namelist.seaice streams.seaice output_%s_%s_%s_%i" %(gridType, operatorMethod, grid, subcycleNumber))
                    os.system("ln -s namelist.seaice.%s.%i namelist.seaice" %(operatorMethod, subcycleNumber))
                    os.system("ln -s streams.seaice.square streams.seaice")

                    nProcs = 1
                    execute_model(MPAS_SEAICE_EXECUTABLE,
                                  nProcs)

                    os.system("mv output output_%s_%s_%s_%i" %(gridType, operatorMethod, grid, subcycleNumber))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
