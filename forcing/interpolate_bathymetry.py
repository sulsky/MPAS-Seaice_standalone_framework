from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
import subprocess
import sys
import glob
import configparser
from create_forcing import create_scrip_grid_file, get_mpas_grid_info, create_scrip_file_MPAS, write_scrip_in_file, get_remapping_data
import argparse

#-------------------------------------------------------------------------------

def create_T62_remap_file(filenameT62,
                          filenameScrip,
                          title,
                          latDimname,
                          lonDimname,
                          latVarname,
                          lonVarname):

    fileT62 = Dataset(filenameT62,"r")

    nLat = len(fileT62.dimensions[latDimname])
    nLon = len(fileT62.dimensions[lonDimname])

    LATin = fileT62.variables[latVarname][:]
    LONin = fileT62.variables[lonVarname][:]

    fileT62.close()

    LAT = np.zeros(nLat+2)
    LON = np.zeros(nLon+2)

    LAT[1:-1] = LATin[:]
    LON[1:-1] = LONin[:]

    # center
    latCenter = np.zeros((nLat,nLon))
    lonCenter = np.zeros((nLat,nLon))

    for iLon in range(0,nLon):
       latCenter[:,iLon] = LATin[:]

    for iLat in range(0,nLat):
       lonCenter[iLat,:] = LONin[:]

    latCenter = np.radians(latCenter)
    lonCenter = np.radians(lonCenter)

    # corners
    latCorner = np.zeros((nLat,nLon,4))
    lonCorner = np.zeros((nLat,nLon,4))

    for iLon in range(0,nLon):

        iLon2 = iLon + 1

        lonCorner[:,iLon,0] = 0.5 * (LON[iLon2] + LON[iLon2-1])
        lonCorner[:,iLon,1] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,2] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,3] = 0.5 * (LON[iLon2] + LON[iLon2-1])

    lonCorner = (np.where(lonCorner < 0.0, lonCorner + 360.0, lonCorner))
    lonCorner = np.radians(lonCorner)

    for iLat in range(0,nLat):

        iLat2 = iLat + 1

        latCorner[iLat,:,0] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,1] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,2] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])
        latCorner[iLat,:,3] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])

    latCorner = np.radians(latCorner)

    # create file
    nGridSize = nLat * nLon
    nGridCorners = 4
    gridRank = 2
    gridDims = np.array([nLon,nLat])
    gridImask = np.ones(nGridSize,dtype="i")

    latCornerScrip = np.zeros((nLat*nLon,4))
    lonCornerScrip = np.zeros((nLat*nLon,4))

    for iLat in range(0,nLat):
        for iLon in range(0,nLon):
            for iCorner in range(0,4):
                ij = iLat * nLon + iLon
                latCornerScrip[ij,iCorner] = latCorner[iLat,iLon,iCorner]
                lonCornerScrip[ij,iCorner] = lonCorner[iLat,iLon,iCorner]

    create_scrip_grid_file(filenameScrip,
                           nGridSize,
                           nGridCorners,
                           gridRank,
                           gridDims,
                           latCenter.flatten(),
                           lonCenter.flatten(),
                           gridImask,
                           latCornerScrip,
                           lonCornerScrip,
                           title)

    return nGridSize

#-------------------------------------------------------------------------------

def create_bathymetry(inputDir,
                      filenameIn,
                      varnameIn,
                      outputDir,
                      filenameOut,
                      varnameOut,
                      remapMatrix,
                      dstGridSize):

    # create output file
    filenameOut = outputDir+"/"+filenameOut
    fileForcing = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    # dimensions
    nCells = fileForcing.createDimension("nCells",dstGridSize)

    print("    Variable: %s to %s" %(varnameIn, varnameOut))

    # open input file
    fileInput = Dataset(filenameIn,"r")
    arrayIn = fileInput.variables[varnameIn][:]
    fileInput.close()

    arrayOut = np.zeros(dstGridSize)

    arrayIn = arrayIn[:,:].flatten()
    arrayOut[:] = remapMatrix.dot(arrayIn)

    # output variable to netcdf file
    var = fileForcing.createVariable(varnameOut,"d",dimensions=["nCells"])
    var[:] = arrayOut[:]

    # close forcing file
    fileForcing.close()

#-------------------------------------------------------------------------------

def write_scrip_in_file(srcTitle):

    scripFile = open("scrip_in","w")

    scripFile.write("&remapInputs\n")
    scripFile.write("    num_maps = 1\n")
    scripFile.write("    gridFile1 = 'remap_grid_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    gridFile2 = 'remap_grid_MPAS_tmp.nc'\n")
    scripFile.write("    interpFile1 = 'remap_%s_to_MPAS_tmp.nc'\n" %(srcTitle))
    scripFile.write("    interpFile2 = 'remap_MPAS_to_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    mapName1 = '%s to MPAS bilinear mapping'\n" %(srcTitle))
    scripFile.write("    mapName2 = 'MPAS to %s bilinear mapping'\n" %(srcTitle))
    scripFile.write("    mapMethod = 'bilinear'\n")
    scripFile.write("    normalizeOpt = 'fracArea'\n")
    scripFile.write("    outputFormat = 'scrip'\n")
    scripFile.write("    restrict_type = 'latitude'\n")
    scripFile.write("    num_srch_bins = 90 \n")
    scripFile.write("    luse_grid1_area = .false.\n")
    scripFile.write("    luse_grid2_area = .false.\n")
    scripFile.write("/\n")

    scripFile.close()

#-------------------------------------------------------------------------------

def perform_remapping(meshDir,
                      meshFilename,
                      latDimname,
                      lonDimname,
                      latVarname,
                      lonVarname,
                      inputDir,
                      filenameIn,
                      varnameIn,
                      filenameMPASGrid,
                      outputDir,
                      filenameOut,
                      varnameOut):

    # create MPAS scrip grid file
    print("create_scrip_file_MPAS")
    scripGridFilename  = "remap_grid_MPAS_tmp.nc"
    dstGridSize = create_scrip_file_MPAS(filenameMPASGrid,
                                         scripGridFilename)

    # create T62 remapping file
    print("create_T62_remap_file")
    filenameT62 = meshDir+meshFilename
    scripT62Filename = "remap_grid_T62_tmp.nc"
    srcGridSize = create_T62_remap_file(filenameT62,
                                        scripT62Filename,
                                        "T62",
                                        latDimname,
                                        lonDimname,
                                        latVarname,
                                        lonVarname)

    # create input scrip file
    print("write_scrip_in_file")
    write_scrip_in_file("T62")

    # run scrip to generate weights
    print("ESMF_RegridWeightGen")
    process = subprocess.Popen(["ESMF_RegridWeightGen",
                                "--source",     "remap_grid_T62_tmp.nc",
                                "--destination","remap_grid_MPAS_tmp.nc",
                                "--weight",     "remap_T62_to_MPAS_tmp.nc",
                                "--method",     "bilinear",
                                "--weight_only"])
    process.wait()
    if (process.returncode != 0):
        print("ESMF_RegridWeightGen error: ", process.returncode)
        exit(1);

    # get remapping weights
    print("get_remapping_data")
    filenameRemapping = "remap_T62_to_MPAS_tmp.nc"
    remapMatrix = get_remapping_data(filenameRemapping,
                                     srcGridSize,
                                     dstGridSize)

    # combined output file
    print("create_bathymetry")
    create_bathymetry(inputDir,
                      filenameIn,
                      varnameIn,
                      outputDir,
                      filenameOut,
                      varnameOut,
                      remapMatrix,
                      dstGridSize)

#-------------------------------------------------------------------------------

def create_bathymetry_from_config(configFilename):

    config = configparser.ConfigParser()
    config.read(configFilename)

    # input_mesh
    meshDir          = config.get   ('input_mesh',   'meshDir')
    meshFilename     = config.get   ('input_mesh',   'meshFilename')
    latDimname       = config.get   ('input_mesh',   'latDimname')
    lonDimname       = config.get   ('input_mesh',   'lonDimname')
    latVarname       = config.get   ('input_mesh',   'latVarname')
    lonVarname       = config.get   ('input_mesh',   'lonVarname')

    # input_fields
    inputDir         = config.get   ('input_fields', 'inputDir')
    filenameIn       = config.get   ('input_fields', 'filenameIn')
    varnameIn        = config.get   ('input_fields', 'varnameIn')

    # output
    filenameMPASGrid = config.get   ('output',       'filenameMPASGrid')
    outputDir        = config.get   ('output',       'outputDir')
    filenameOut      = config.get   ('output',       'filenameOut')
    varnameOut       = config.get   ('output',       'varnameOut')

    perform_remapping(meshDir,
                      meshFilename,
                      latDimname,
                      lonDimname,
                      latVarname,
                      lonVarname,
                      inputDir,
                      filenameIn,
                      varnameIn,
                      filenameMPASGrid,
                      outputDir,
                      filenameOut,
                      varnameOut)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create bathymetry file from lat/lon input data')

    parser.add_argument('-c', dest='configFilename', required=True, help='Config filename')

    args = parser.parse_args()

    create_bathymetry_from_config(args.configFilename)
