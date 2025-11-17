import glob
from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def number_of_icebergs():

    filenames = sorted(glob.glob("./output/icebergs_out*"))

    for filename in filenames:

        filein = Dataset(filename,"r")

        try:
            statusIB = filein.variables["statusIB"][0,:]
            nIcebergs = np.sum(statusIB)
        except:
            nIcebergs = 0

        print(filename, nIcebergs)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    number_of_icebergs()
