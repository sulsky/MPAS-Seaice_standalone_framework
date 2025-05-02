import sys

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
particleFilename2 = "/particles_output.2000-01-11_00.00.00.nc"

outDirs = [
    "output_cosine_bell_10242",
    "output_cosine_bell_163842",
    "output_cosine_bell_2562",
    "output_cosine_bell_40962",
    "output_slotted_cylinder_10242",
    "output_slotted_cylinder_163842",
    "output_slotted_cylinder_2562",
    "output_slotted_cylinder_40962"]

usePolympos = [False, True]

print("Get testcase data")
get_testcase_data()

print("Create ICs")
create_ics()

print("Add deldyn to ICs")
add_deldyn_to_ics(3600.0)

print("Create particles")
create_particles()

for usePolympo in usePolympos:

    print("usePolympo: ", usePolympo)

    nmlChanges = {"mpm":{"config_use_mpm_polympo":usePolympo}}
    create_new_namelist("namelist.seaice.advection_convergence", "namelist.seaice", nmlChanges)

    print("Run models")
    run_model(usePolympo)

    print("Check particles moved")
    for outDir in outDirs:
        check_particles_moved(outDir+particleFilename1,
                              outDir+particleFilename2)

    print("Plot test case")
    if (usePolympo):
        runtype = "polympo"
    else:
        runtype = "original"
    plot_testcase(runtype)

    print("Advection map")
    advection_map()

    print("Advection equatorial")
    advection_equatorial()

    print("Advection error convergence")
    advection_error_convergence()
