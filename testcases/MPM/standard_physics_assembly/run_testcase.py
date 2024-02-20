import sys
from run_model import run_model
from check_results import check_results
sys.path.append("../advection")
from check_particle_positions_nprocs import check_particle_positions_nprocs
sys.path.append("../../../testing")
from compare_mpas_files import compare_files
from testing_utils import get_domain, print_colour, create_test_directory, test_summary
import os
sys.path.append("../../../utils/MPM/particle_initialization/")
from empty_particle_file import empty_particle_file

#-------------------------------------------------------------------------------

def check_run(n1, n2):

        # make a test directory
        testDir = "testDir%i_%i" %(n1, n2)
        create_test_directory(testDir)
        os.chdir("../")

        # make log file
        logfile = open("log_test.txt", "w")
        title = "Test: Parallelism, n1=%i, n2=%i" %(n1, n2)
        print_colour(title, "title")
        logfile.write(title)

        # run comparison
        file1="output_%i/output.2000.nc" %(n1)
        file2="output_%i/output.2000.nc" %(n2)
        nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile)
        failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "standard_physics_assembly")

        if (os.path.isfile("vars_differ.nc")):
                cmd = "mv vars_differ.nc %s" %(testDir)
                os.system(cmd)

        check_particle_positions_nprocs(n1,n2)

#-------------------------------------------------------------------------------

# domains directory
domainsDir = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
if (domainsDir == None):
        raise Exception("Environment variable MPAS_SEAICE_DOMAINS_DIR must be set if no domains directory specified")
if (not os.path.exists(domainsDir)):
        raise Exception("Requested domains directory does not exist")


# get domain
domain="domain_QU120km"
get_domain(domainsDir, domain)

# empty particles
empty_particle_file("particles.nc")


# run models
run_model(1)
run_model(16)
run_model(32)


# check runs output
check_run(1,16)
check_run(16,32)
