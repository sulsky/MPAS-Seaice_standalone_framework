import glob
from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def number_of_icebergs():

    filenames = sorted(glob.glob("./output/icebergs_out*"))

    for filename in filenames:

        filein = Dataset(filename,"r")

        nIcebergs = len(filein.dimensions["nIcebergs"])

        print(filename, nIcebergs)


#-------------------------------------------------------------------------------

if __name__ == "__main__":

    number_of_icebergs()
