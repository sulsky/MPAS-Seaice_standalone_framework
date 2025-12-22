import glob
from netCDF4 import Dataset
import numpy as np
import argparse

#-------------------------------------------------------------------------------

def number_of_icebergs(filenameTemplate):

    filenames = sorted(glob.glob(filenameTemplate))

    print("nIcebergsStatus, nIcebergs, nIcebergsCell, nIcebergsPhysical")

    for filename in filenames:

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        try:
            statusIB = filein.variables["statusIB"][0,:]
            nIcebergsStatus = np.sum(statusIB)
        except:
            nIcebergsStatus = 0

        nIcebergsCell = filein.variables["nIcebergsCell"][:]
        nIcebergsCell = np.sum(nIcebergsCell)

        try:
            nIcebergCategories = len(filein.dimensions["nIcebergCategories"])
            icebergCategoryFluxScaling = filein.variables["icebergCategoryFluxScaling"][:]
            icebergCategory = filein.variables["icebergCategory"][0,:]
            nIcebergsPhysical = 0
            for iIceberg in range(0,nIcebergs):
                if (statusIB[iIceberg] == 1):
                    nIcebergsPhysical += icebergCategoryFluxScaling[icebergCategory[iIceberg]-1]
            nIcebergsPhysical = int(nIcebergsPhysical)
        except:
            nIcebergsPhysical = -1

        print(filename, nIcebergsStatus, nIcebergs, nIcebergsCell, nIcebergsPhysical)

        filein.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-f', dest='filenameTemplate', default="./output/icebergs_out*", help='')

    args = parser.parse_args()

    number_of_icebergs(args.filenameTemplate)
