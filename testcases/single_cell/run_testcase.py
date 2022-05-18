import os
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_testcase():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""
    MPAS_SEAICE_DOMAINS_DIR = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')

    # copy namelist and streams file
    os.system("cp ../../configurations/standard_physics_single_cell/namelist.seaice .")
    os.system("cp ../../configurations/standard_physics_single_cell/streams.seaice .")

    # forcing
    os.system("python %s/domain_sc_71.35_-156.5/get_domain.py" %(MPAS_SEAICE_DOMAINS_DIR))

    # run MPAS-Seaice
    os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

    # plot output
    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
