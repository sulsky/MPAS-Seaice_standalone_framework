import subprocess
import os

#-------------------------------------------------------------------------------

def get_testcase_data():

    nCellsArray = [10242, 163842, 2562, 40962]

    dirName = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testcases/strain_stress_divergence/"

    for nCells in nCellsArray:

        # grid file
        filename = "grid.%i.nc" %(nCells)

        if (not os.path.isfile(filename)):

            args = ["wget", dirName+filename]

            process = subprocess.Popen(args, stdout=subprocess.PIPE)

            while process.poll() is None:
                line = process.stdout.readline()

        # graph file
        filename = "graph.%i.info" %(nCells)

        if (not os.path.isfile(filename)):

            args = ["wget", dirName+filename]

            process = subprocess.Popen(args, stdout=subprocess.PIPE)

            while process.poll() is None:
                line = process.stdout.readline()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data()
