import glob
from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def number_of_icebergs(filenameTemplate):

    filenames = sorted(glob.glob(filenameTemplate))

    for filename in filenames:

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        try:
            statusIB = filein.variables["statusIB"][0,:]
            nIcebergsStatus = np.sum(statusIB)
        except:
            nIcebergsStatus = 0

        print(filename, nIcebergsStatus, nIcebergs)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', default="./output/icebergs_out*", help='')

    args = parser.parse_args()

    number_of_icebergs(args.filenameTemplate)
