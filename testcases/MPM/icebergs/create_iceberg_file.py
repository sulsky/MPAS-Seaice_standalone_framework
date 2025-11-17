from netCDF4 import Dataset
import numpy as np
from math import radians

#-------------------------------------------------------------------------------

def write_iceberg_file(nIcebergs,
                       latIceberg,
                       lonIceberg,
                       icebergLength,
                       icebergHeight):

    fileout = Dataset("icebergs.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nIcebergs",nIcebergs)

    var = fileout.createVariable("latIceberg", "d", dimensions=["nIcebergs"])
    var[:] = latIceberg[:]

    var = fileout.createVariable("lonIceberg", "d", dimensions=["nIcebergs"])
    var[:] = lonIceberg[:]

    var = fileout.createVariable("icebergLength", "d", dimensions=["nIcebergs"])
    var[:] = icebergLength[:]

    var = fileout.createVariable("icebergHeight", "d", dimensions=["nIcebergs"])
    var[:] = icebergHeight[:]

    fileout.close()

#-------------------------------------------------------------------------------

def create_iceberg_file():

    nIcebergs = 1

    latIceberg = np.array([radians(-72.216494)])
    lonIceberg = np.array([radians(-19.975072)])
    icebergLength = np.array([2000.0])
    icebergHeight = np.array([73.0])

    write_iceberg_file(nIcebergs,
                       latIceberg,
                       lonIceberg,
                       icebergLength,
                       icebergHeight)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_iceberg_file()
