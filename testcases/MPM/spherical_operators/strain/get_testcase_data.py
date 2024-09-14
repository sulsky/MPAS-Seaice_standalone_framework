import subprocess
from os.path import exists

#-------------------------------------------------------------------------------

def get_testcase_data():

    filenames = ["grid.10242.nc",
                 "grid.163842.nc",
                 "grid.2562.nc",
                 "grid.40962.nc"]

    dirName = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testcases/strain_stress_divergence/"

    for filename in filenames:

        totalFilename = dirName+filename

        if (not exists(filename)):

            args = ["wget", totalFilename]

            process = subprocess.Popen(args, stdout=subprocess.PIPE)

            while process.poll() is None:
                line = process.stdout.readline()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data()
