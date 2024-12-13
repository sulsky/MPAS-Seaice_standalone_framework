from create_grids import create_grids
from create_ics import create_ics
from create_particles import create_particles
from run_model1 import run_model
from average_variational_strains import average_variational_strains
from average_weak_strains import average_weak_strains
from strain_map import strain_map
from strain_scaling1 import strain_scaling
from velocity_scaling import velocity_scaling
from error_analysis_strain import error_analysis_strain

create_grids()

create_ics()

create_particles()

run_model()

#average_variational_strains()

#average_weak_strains()

#strain_map()

strain_scaling()

error_analysis_strain()

velocity_scaling()
