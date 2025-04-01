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
particleFilename2 = "/particles_output.2000-03-01_00.00.00.nc"

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

get_testcase_data()

create_ics()

add_deldyn_to_ics(3600.0)

create_particles()

for usePolympo in usePolympos:

    print("usePolympo: ", usePolympo)

    nmlChanges = {"mpm":{"config_use_mpm_polympo":usePolympo}}
    create_new_namelist("namelist.seaice.advection_convergence", "namelist.seaice", nmlChanges)

    run_model(usePolympo)

    for outDir in outDirs:
        check_particles_moved(outDir+particleFilename1,
                              outDir+particleFilename2)

    if (usePolympo):
        runtype = "polympo"
    else:
        runtype = "original"
    plot_testcase(runtype)

    advection_map()

    advection_equatorial()

    advection_error_convergence()
