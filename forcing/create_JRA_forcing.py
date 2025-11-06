import os
import uuid
import argparse
import calendar
from create_atmos_forcing import create_forcing_from_config

varnames = ["ugrd10m","vgrd10m","tmp2m","spfh2m"]

filenamePrefixes = {"ugrd10m":"jra3q.anl_surf125.0_2_2.ugrd10m-hgt-an-ll125",
                    "vgrd10m":"jra3q.anl_surf125.0_2_3.vgrd10m-hgt-an-ll125",
                    "tmp2m"  :"jra3q.anl_surf125.0_0_0.tmp2m-hgt-an-ll125",
                    "spfh2m" :"jra3q.anl_surf125.0_1_0.spfh2m-hgt-an-ll125"}

daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]

jraURL = "https://data.rda.ucar.edu/d640000/"

#-------------------------------------------------------------------------------

def download_input_data(tmpDirname,
                        year):

    for varname in varnames:
        for month in range(0,12):

            daysInThisMonth = daysInMonth[month]
            if (calendar.isleap(year) and month == 1):
                daysInThisMonth += 1

            directory = "anl_surf125/%4.4i%2.2i/" %(year,month+1)
            filenamePrefix = filenamePrefixes[varname]
            dateStart = "%4.4i%2.2i0100" %(year,month+1)
            dateEnd   = "%4.4i%2.2i%2.2i18" %(year,month+1,daysInThisMonth)
            filename = "%s.%s_%s.nc" %(filenamePrefix,dateStart,dateEnd)

            cmd = "wget -P %s%s %s%s%s" %(tmpDirname,varname,jraURL,directory,filename)
            print(cmd)
            os.system(cmd)

#-------------------------------------------------------------------------------

def concatenate_input_data(tmpDirname,
                           year):

    for varname in varnames:

        os.chdir("%s/%s" %(tmpDirname,varname))

        outputFilename = "%s.%i.nc" %(varname,year)
        cmd = "ncrcat * -o %s" %(outputFilename)
        print(cmd)
        os.system(cmd)

        os.chdir("../../")

#-------------------------------------------------------------------------------

def create_config_file(filenameConfig,
                       year,
                       inputDir,
                       mpasMeshFilename,
                       outputDir):

    fileout = open(filenameConfig,"w")

    fileout.write("[time]\n")
    fileout.write("startYear = %i\n" %(year))
    fileout.write("endYear = %i\n" %(year))
    fileout.write("timeType = sixhourly_gregorian\n")
    fileout.write("\n")
    fileout.write("[input_mesh]\n")
    fileout.write("meshDir = %s\n" %(inputDir))
    fileout.write("meshFilename = /tmp2m/tmp2m.%i.nc\n" %(year))
    fileout.write("latDimname = lat\n")
    fileout.write("lonDimname = lon\n")
    fileout.write("latVarname = lat\n")
    fileout.write("lonVarname = lon\n")
    fileout.write("\n")
    fileout.write("[input_fields]\n")
    fileout.write("inputDir = %s\n" %(inputDir))
    fileout.write("filenameTemplatesIn =\n")
    fileout.write("                    /tmp2m/tmp2m.$Y.nc\n")
    fileout.write("                    /spfh2m/spfh2m.$Y.nc\n")
    fileout.write("                    /ugrd10m/ugrd10m.$Y.nc\n")
    fileout.write("                    /vgrd10m/vgrd10m.$Y.nc\n")
    fileout.write("varnamesIn =\n")
    fileout.write("           tmp2m-hgt-an-ll125\n")
    fileout.write("           spfh2m-hgt-an-ll125\n")
    fileout.write("           ugrd10m-hgt-an-ll125\n")
    fileout.write("           vgrd10m-hgt-an-ll125\n")
    fileout.write("\n")
    fileout.write("[output]\n")
    fileout.write("filenameMPASGrid = %s\n" %(mpasMeshFilename))
    fileout.write("outputDir = %s\n" %(outputDir))
    fileout.write("filenameOutTemplate = LYq_six_hourly.$Y.nc\n")
    fileout.write("varnamesOut =\n")
    fileout.write("            airTemperature\n")
    fileout.write("            airSpecificHumidity\n")
    fileout.write("            uAirVelocity\n")
    fileout.write("            vAirVelocity\n")

    fileout.close()

#-------------------------------------------------------------------------------

def create_JRA_forcing(mpasMeshFilename,
                       yearStart,
                       yearEnd,
                       outputDir):

    tmpDirname = "tmp_%s/" %(uuid.uuid4().hex)
    os.mkdir(tmpDirname)

    for year in range(yearStart, yearEnd+1):

        # download input data
        download_input_data(tmpDirname,
                            year)

        # concatenate input data
        concatenate_input_data(tmpDirname,
                               year)

        # create forcing config file
        filenameConfig = tmpDirname + "config_jra_forcing.cfg"
        create_config_file(filenameConfig,
                           year,
                           tmpDirname,
                           mpasMeshFilename,
                           outputDir)

        # create forcing file
        create_forcing_from_config(filenameConfig)

        # clean up
        cmd = "rm -rf %s*" %(tmpDirname)
        print(cmd)
        if (len(tmpDirname) == 37):
            os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m',   dest='mpasMeshFilename', required=True, help='')
    parser.add_argument('--y0', dest='yearStart',        required=True, type=int, help='')
    parser.add_argument('--y1', dest='yearEnd',          required=True, type=int, help='')
    parser.add_argument('-o',   dest='outputDir',        required=True, help='')

    args = parser.parse_args()

    create_JRA_forcing(args.mpasMeshFilename,
                       args.yearStart,
                       args.yearEnd,
                       args.outputDir)
