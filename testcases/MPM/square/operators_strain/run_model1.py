import os
try:
    import f90nml
except ImportError:
    print("Module f90nml needed and not available")
    raise

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    operatorMethods = ["mpmvar","mpmweak"]

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
    #grids = {"hex" :["0656x0752"]}

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

                if (operatorMethod == "mpmvar"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                    "config_stress_divergence_scheme":"variational"},
                                "use_sections": {"config_use_mpm": True}}
                elif (operatorMethod == "mpmweak"):
                    nmlPatch = {"velocity_solver": {"config_strain_scheme":"mpm",
                                                    "config_stress_divergence_scheme":"weak"},
                                "use_sections": {"config_use_mpm": True}}

                f90nml.patch("namelist.seaice.strain", nmlPatch, "namelist.seaice.%s" %(operatorMethod))

                os.system("rm -rf namelist.seaice streams.seaice output_%s_%s_%s" %(gridType, operatorMethod, grid))
                os.system("ln -s namelist.seaice.%s namelist.seaice" %(operatorMethod))
                os.system("ln -s streams.seaice.strain streams.seaice")

                os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

                os.system("mv output output_%s_%s_%s" %(gridType, operatorMethod, grid))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
