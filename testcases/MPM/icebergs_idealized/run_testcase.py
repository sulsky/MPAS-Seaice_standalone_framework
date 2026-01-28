import sys

sys.path.append("../../../utils/testcases")
from get_testcase_data_spherical import get_testcase_data_spherical

from plot_results import plot_testcase

import numpy as np
from netCDF4 import Dataset
from math import pi
import os

#-------------------------------------------------------------------------------

def write_iceberg_file(nIcebergs,
                       latIceberg,
                       lonIceberg,
                       icebergLength,
                       icebergHeight,
                       filename):

    fileout = Dataset(filename,"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nIcebergs",nIcebergs)

    var = fileout.createVariable("cellIDCreationIB", "i", dimensions=["nIcebergs"])
    var[:] = 1

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

def write_ic_file_air_drag(nCells):

    fileout = Dataset("ic_air_drag.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 1.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_ocean_drag(nCells):

    fileout = Dataset("ic_ocean_drag.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 1.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_coriolis():

    fileout = Dataset("ic_coriolis.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nIcebergs",1)

    var = fileout.createVariable("uVelocityIceberg","d",dimensions=["nIcebergs"])
    var[:] = 0.0

    var = fileout.createVariable("vVelocityIceberg","d",dimensions=["nIcebergs"])
    var[:] = 1.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_surface_tilt(nCells):

    fileout = Dataset("ic_surface_tilt.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 1.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_wave_radiation(nCells):

    fileout = Dataset("ic_wave_radiation.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)

    var = fileout.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 1.0

    var = fileout.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 0.0

    var = fileout.createVariable("uOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = -1.0

    var = fileout.createVariable("vOceanVelocityIcebergCell","d",dimensions=["nCells"])
    var[:] = 0.0

    fileout.close()

#-------------------------------------------------------------------------------

def write_ic_file_seaice_force(nCells,
                               nVertices,
                               concLabel,
                               iceConcentration,
                               icePressure):

    fileout = Dataset("ic_seaice_force_%s.nc" %(concLabel),"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nVertices",nVertices)

    var = fileout.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = 1.0

    var = fileout.createVariable("vVelocity","d",dimensions=["nVertices"])
    var[:] = 0.0

    var = fileout.createVariable("iceAreaCell","d",dimensions=["nCells"])
    var[:] = iceConcentration

    var = fileout.createVariable("icePressure","d",dimensions=["nCells"])
    var[:] = icePressure

    fileout.close()

#-------------------------------------------------------------------------------

def run_testcase():

    nCells = 2562

    print("Get testcase data")
    print("=================")
    get_testcase_data_spherical()

    print("Grid data")
    print("=========")
    try:
        os.remove("grid.nc")
    except:
        pass
    os.symlink("grid.%i.nc" %(nCells), "grid.nc")

    filenameGrid = "grid.nc"
    fileGrid = Dataset(filenameGrid,"r")
    nVertices = len(fileGrid.dimensions["nVertices"])
    fileGrid.close()

    print("Iceberg particle input")
    print("======================")
    write_iceberg_file(1,
                       np.array([0.0]),
                       np.array([0.0]),
                       np.array([1000.0]),
                       np.array([1000.0]),
                       "icebergs_equator.nc")

    write_iceberg_file(1,
                       np.array([pi/4.0]),
                       np.array([0.0]),
                       np.array([1000.0]),
                       np.array([1000.0]),
                       "icebergs_midlatitude.nc")

    print("Iceberg IC data")
    print("===============")
    write_ic_file_air_drag(nCells)
    write_ic_file_ocean_drag(nCells)
    write_ic_file_surface_tilt(nCells)
    write_ic_file_wave_radiation(nCells)
    write_ic_file_seaice_force(nCells, nVertices, "low_conc",  0.01, 0.0)
    write_ic_file_seaice_force(nCells, nVertices, "mid_conc",  0.5,  0.0)
    write_ic_file_seaice_force(nCells, nVertices, "high_conc", 0.99, 1.0e5)
    write_ic_file_coriolis()

    testcases = ["air_drag",
                 "ocean_drag",
                 "coriolis",
                 "surface_tilt",
                 "wave_radiation",
                 "seaice_force_low_conc",
                 "seaice_force_mid_conc",
                 "seaice_force_high_conc"]

    print("Individual test cases")
    print("=====================")
    for testcase in testcases:

        print("   ", testcase)

        try:
            os.remove("namelist.seaice")
        except:
            pass
        try:
            os.remove("streams.seaice")
        except:
            pass

        os.symlink("namelist.seaice.%s" %(testcase),"namelist.seaice")
        os.symlink("streams.seaice.%s" %(testcase),"streams.seaice")

        cmd = "../../../../MPAS-Seaice-MPM/components/mpas-seaice/seaice_model"
        os.system(cmd)

        os.rename("output", "output_%s" %(testcase))

        plot_testcase(testcase)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
