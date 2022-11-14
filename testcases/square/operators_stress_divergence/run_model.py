import os
import argparse
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model(testName, ignoreWeak=False):

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    if (not ignoreWeak):
        operatorMethods = ["wachspress","pwl","weak"]
    else:
        operatorMethods = ["wachspress","pwl"]

    gridTypes = ["hex","quad"]
    #gridTypes = ["hex"]

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

    for gridType in gridTypes:

        print("Grid type: ", gridType)

        for operatorMethod in operatorMethods:

            print("  Operator Method: ", operatorMethod)

            for grid in grids[gridType]:

                print("    Grid: ", grid)

                os.system("rm grid.nc")
                os.system("rm ic.nc")

                os.system("ln -s grid_%s_%s_%s.nc grid.nc" %(testName,gridType,grid))
                os.system("ln -s ic_%s_%s_%s.nc ic.nc" %(testName,gridType,grid))

                if (operatorMethod == "wachspress"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                    "config_stress_divergence_scheme":"variational",
                                                    "config_variational_basis":"wachspress",
                                                    "config_variational_denominator_type":"alternate"}}
                elif (operatorMethod == "pwl"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"variational",
                                                    "config_stress_divergence_scheme":"variational",
                                                    "config_variational_basis":"pwl",
                                                    "config_variational_denominator_type":"alternate"}}
                elif (operatorMethod == "weak"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"weak",
                                                    "config_stress_divergence_scheme":"weak"}}

                f90nml.patch("namelist.seaice.stress_divergence", nmlPatch, "namelist.seaice.%s" %(operatorMethod))

                os.system("rm -rf namelist.seaice streams.seaice output_%s_%s_%s" %(gridType, operatorMethod, grid))
                os.system("ln -s namelist.seaice.%s namelist.seaice" %(operatorMethod))
                os.system("ln -s streams.seaice.stress_divergence streams.seaice")

                os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

                outputDir = "output_%s_%s_%s_%s" %(testName, gridType, operatorMethod, grid)
                os.system("mv output %s" %(outputDir))
                os.system("cp log.seaice.0000.out %s" %(outputDir))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-t', dest='testName', help='')
    parser.add_argument('-w', dest='ignoreWeak', action='store_true', help='')

    args = parser.parse_args()

    run_model(args.testName, args.ignoreWeak)
