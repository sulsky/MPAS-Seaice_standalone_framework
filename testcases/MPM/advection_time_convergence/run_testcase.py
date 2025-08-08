import os, sys, math, f90nml, argparse

from run_model import run_model
from plot_testcase import plot_testcase
from advection_map import advection_map
from advection_equatorial import advection_equatorial
from advection_error_convergence import advection_error_convergence
from create_particles import create_particles

sys.path.append("../../advection")
from create_ics import create_ic_file

sys.path.append("../advection")
from add_deldyn_to_ics import add_deldyn_to_ics
from check_particles_moved import check_particles_moved

sys.path.append("../../../testing")
from testing_utils import create_new_namelist

sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

#-------------------------------------------------------------------------------

def run_testcase(logFilename=None):

    if (logFilename is not None):
        logFile = open(logFilename,"a")
        logFile.write("\nAdvection test case\n")
        logFile.write(  "===================\n")
        logFile.flush()
    else:
        logFile = None

    particleFilename1 = "/particles_output.2000-01-01_00.00.00.nc"
    particleFilename2 = "/particles_output.2000-01-06_00.00.00.nc"

    outDirs = [
        "output_cosine_bell_3600",
        "output_cosine_bell_1800",
        "output_cosine_bell_900",
        "output_cosine_bell_450",
        "output_slotted_cylinder_3600",
        "output_slotted_cylinder_1800",
        "output_slotted_cylinder_900",
        "output_slotted_cylinder_450"]

    res = "2562"
    icTypes = ["cosine_bell", "slotted_cylinder"]
    DynamicsTimeStep = [3600, 1800, 900, 450]
    usePolympos = [False, True]

    print("Get testcase data")
    print("=================")
    get_testcase_data_spherical()

    print("\nCreate ICs")
    print("==========")
    create_ic_file(res, "slotted_cylinder" , math.pi / 6.0)
    create_ic_file(res, "cosine_bell"      , math.pi / 6.0)

    print("\nAdd deldyn to IC")
    print("==========")
    for icType in icTypes:
       cmd = "cp ic_%s_%s.nc ic_%s.nc" %(icType,res,icType)
       os.system(cmd)
    for tstep in DynamicsTimeStep:
       add_deldyn_to_ics(tstep,res)
       for icType in icTypes:
            cmd = "cp ic_%s_%s.nc ic_%s_%i.nc" %(icType,res,icType,tstep)
            os.system(cmd)
            cmd = "cp ic_%s.nc ic_%s_%s.nc" %(icType,icType,res)
            os.system(cmd)

    print("\nCreate particles")
    print("================")
    create_particles()

    print("\nCreate namelists")
    print("================")

    for usePolympo in usePolympos:

        print("usePolympo: ", usePolympo)
        if (usePolympo):
             usePolympoStr = "polympo"
        else:
             usePolympoStr = "nonpolympo"

        t = 0
        for tstep in DynamicsTimeStep:

            nmlChanges = {"mpm":{"config_use_mpm_polympo":usePolympo},
                          "seaice_model":{"config_dt":DynamicsTimeStep[t]}}

            new_namelist = "namelist.seaice.%s.%i" %(usePolympoStr, tstep)
            create_new_namelist("namelist.seaice.advection", new_namelist, nmlChanges)

            t = t + 1

        print("\nRun models")
        print("==========")
        run_model(res,
                  1,
                  usePolympoStr,
                  logFile)

        print("\nCheck particles moved")
        print("=====================")
        for outDir in outDirs:
            check_particles_moved(outDir+"_"+usePolympoStr+particleFilename1,
                                  outDir+"_"+usePolympoStr+particleFilename2)

        print("\nPlot test case")
        print("==============")
        if (usePolympo):
            runtype = "polympo"
        else:
            runtype = "nonpolympo"
        plot_testcase(runtype)

        print("\nAdvection map")
        print("=============")
        advection_map(runtype)

        print("\nAdvection equatorial")
        print("====================")
        advection_equatorial(runtype)

        print("\nAdvection error convergence")
        print("===========================")
        advection_error_convergence(runtype)

    if (logFile is not None):
        logFile.close()


#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-l', dest='logFilename')

    args = parser.parse_args()

    run_testcase(args.logFilename)
