from create_particles import create_particles
from run_model import run_model
from average_variational_strains import average_variational_strains
from average_weak_strains import average_weak_strains
from strain_map import strain_map
from strain_scaling import strain_scaling
from velocity_scaling import velocity_scaling
from error_analysis_strain import error_analysis_strain
import sys
sys.path.append("../../../square/operators_strain")
from create_grids import create_grids
from create_ics import create_ics

create_grids()

create_ics()

create_particles()

run_model()

strain_scaling()

error_analysis_strain()

velocity_scaling()
