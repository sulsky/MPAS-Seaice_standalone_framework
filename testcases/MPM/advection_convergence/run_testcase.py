import sys, math, f90nml

sys.path.append("../../advection")

from get_testcase_data import get_testcase_data
from create_ics import create_ics
from create_particles import create_particles
from run_model import run_model
from plot_testcase import plot_testcase
from advection_map import advection_map
from advection_equatorial import advection_equatorial
from advection_error_convergence import advection_error_convergence

sys.path.append("../advection")
from add_deldyn_to_ics import add_deldyn_to_ics
from check_particles_moved import check_particles_moved

sys.path.append("../../../testing")
from testing_utils import create_new_namelist

particleFilename1 = "/particles_output.2000-01-01_00.00.00.nc"
particleFilename2 = "/particles_output.2000-01-06_00.00.00.nc"

outDirs = [
    "output_cosine_bell_10242",
    "output_cosine_bell_163842",
    "output_cosine_bell_2562",
    "output_cosine_bell_40962",
    "output_slotted_cylinder_10242",
    "output_slotted_cylinder_163842",
    "output_slotted_cylinder_2562",
    "output_slotted_cylinder_40962"]

reses = ["2562", "10242", "40962", "163842"]
icTypes = ["cosine_bell", "slotted_cylinder"]
DynamicsTimeStep = [3600.0, 1800.0, 900.0, 450.0]
usePolympos = [False, True]

print("Get testcase data")
print("=================")
get_testcase_data()

print("\nCreate ICs")
print("==========")
create_ics(math.pi / 6.0)

print("\nAdd deldyn to IC")
print("==========")
for icType in icTypes:
    r = 0
    for res in reses:
        add_deldyn_to_ics(DynamicsTimeStep[r],res)
        r = r + 1

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

    r = 0
    for res in reses:

        nmlChanges = {"mpm":{"config_use_mpm_polympo":usePolympo},
                      "seaice_model":{"config_dt":DynamicsTimeStep[r]}}

        new_namelist = "namelist.seaice.%s.%s" %(usePolympoStr, res)
        create_new_namelist("namelist.seaice.advection_convergence", new_namelist, nmlChanges)

        r = r + 1

    print("\nRun models")
    print("==========")
    run_model(usePolympo)

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
