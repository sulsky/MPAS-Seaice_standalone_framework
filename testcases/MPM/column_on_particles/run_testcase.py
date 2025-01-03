import os
import sys
sys.path.append("../../../testing")
from testing_utils import get_domain, add_pio_namelist_changes, create_new_namelist, test_summary
sys.path.append("../../../utils/MPM/particle_initialization/")
from empty_particle_file import empty_particle_file
from convert_particles_to_cell import convert_particles_to_cell
from compare_mpas_files import compare_files
import argparse

#-------------------------------------------------------------------------------

def add_output_fields_to_stream(filenameIn,
                                filenameOut,
                                changes,
                                fieldsToAdd):

    try:
        import xml.etree.ElementTree as ET
    except ImportError:
        print("Module xml.etree.ElementTree needed and not available")
        sys.exit()

    tree = ET.parse(filenameIn)
    streams = tree.getroot()

    for stream in streams:

        for change in changes:
            if (stream.get('name') == change["streamName"]):
                stream.set(change["attributeName"],change["newValue"])

        for streamsToModify in fieldsToAdd:
            streamName = streamsToModify["streamName"]
            fieldNames = streamsToModify["fieldNames"]
            if (stream.get('name') == streamName):

                for field in list(stream):
                    stream.remove(field)

                for fieldName in fieldNames:
                    field = ET.SubElement(stream, "var")
                    field.set("name",fieldName)

    tree.write(filenameOut)

#-------------------------------------------------------------------------------

def run_testcase(nProcs,
                 runDuration,
                 outputInterval):

    fieldNamesCell = [
        "iceAreaCell",
        "iceVolumeCell",
        "snowVolumeCell",
        "iceAreaCategory",
        "iceVolumeCategory",
        "snowVolumeCategory",
        "surfaceTemperatureCell",
        "shortwaveDown",
        "longwaveDown",
        "seaSurfaceSalinity",
        "shortwaveScalingFactor",
        "airTemperature",
        "congelation",
        "frazilFormation",
        "snowiceFormation",
        "snowMelt",
        "surfaceIceMelt",
        "basalIceMelt",
        "lateralIceMelt"]

    fieldNamesParticle = []
    for fieldName in fieldNamesCell:
        fieldNamesParticle.append(fieldName+"MP")

    fieldNamesCell.append("xtime")
    fieldNamesParticle.append("xtime")

    logfile = open("log_test.txt","w")
    logfile.write("Column physics on particles test\n")

    # get domain
    logfile.write("Set up domains\n")
    domainsDir = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
    if (domainsDir == None):
        raise Exception("Environment variable MPAS_SEAICE_DOMAINS_DIR must be set if no domains directory specified")
    if (not os.path.exists(domainsDir)):
        raise Exception("Requested domains directory does not exist")

    domain = "domain_QU120km"
    get_domain(domainsDir, domain)

    # empty particles
    empty_particle_file("particles.nc")

    # executable
    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")

    # create output file if doesn't exist
    if (not os.path.isdir("output")):
        os.mkdir("output")

    # column in cells run
    logfile.write("Column on cells run\n")

    cmd = "cp ../../../configurations/standard_physics/namelist.seaice namelist.seaice.original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    logfile.write("  Modify namelist\n")
    nmlChanges = {"seaice_model":{"config_run_duration":runDuration}}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)
    create_new_namelist("namelist.seaice.original", "namelist.seaice", nmlChanges)

    cmd = "cp ../../../configurations/standard_physics/streams.seaice streams.seaice.original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    changes = [{"streamName":'output',
                "attributeName":'output_interval',
                "newValue":outputInterval}]
    fieldsToAdd = [{"streamName":'output',
                    "fieldNames":fieldNamesCell}]
    add_output_fields_to_stream("streams.seaice.original",
                                "streams.seaice",
                                changes,
                                fieldsToAdd)

    cmd = "mpirun -np %i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE)
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    cmd = "rm -rf output_original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    cmd = "mv output output_original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    cmd = "mv log.seaice.0000.out log.seaice.0000.out_original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    # column on particles run
    logfile.write("Column on particles run\n")

    cmd = "cp ../../../configurations/standard_physics/namelist.seaice namelist.seaice.original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    logfile.write("  Modify namelist\n")
    nmlChanges = {"seaice_model":{"config_run_duration":runDuration},
                  "mpm": {"config_use_mpm_velocity_solver":False,
                          "config_use_mpm_transport":False,
                          "config_mpm_particle_init_type":'cell',
                          "config_mpm_particle_posn_init_type":'cell'},
                  "use_sections": {"config_use_mpm":True},
                  "column_package": {"config_column_element_type":'particles'}}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)
    create_new_namelist("namelist.seaice.original", "namelist.seaice", nmlChanges)

    cmd = "cp ../../../configurations/standard_physics/streams.seaice namelist.seaice.original"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    changes = [{"streamName":'output',
                "attributeName":'output_interval',
                "newValue":outputInterval}]
    fieldsToAdd = [{"streamName":'output',
                    "fieldNames":fieldNamesParticle}]
    add_output_fields_to_stream("streams.seaice.original",
                                "streams.seaice",
                                changes,
                                fieldsToAdd)

    cmd = "mpirun -np %i %s" %(nProcs, MPAS_SEAICE_EXECUTABLE)
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    cmd = "rm -rf output_mpm"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    cmd = "mv output output_mpm"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    cmd = "mv log.seaice.0000.out log.seaice.0000.out_mpm"
    logfile.write("  %s\n" %(cmd))
    print(cmd)
    os.system(cmd)

    # convert particle column data to cells
    logfile.write("Convert particle files to cell format for comparison\n")
    convert_particles_to_cell("./output_mpm/output.2000.nc",
                              "./output_mpm/output.2000.particles.nc")

    # comparison
    logfile.write("Output file comparison\n")
    file1 = "./output_original/output.2000.nc"
    file2 = "./output_mpm/output.2000.particles.nc"
    logfile.write("  file1: %s\n" %(file1))
    logfile.write("  file2: %s\n" %(file2))

    cmd = "rm vars_differ.nc"
    os.system(cmd)
    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "particles_on_column")

    logfile.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', dest="nProcs", type=int, default=1)
    parser.add_argument('-d', dest="runDuration", default='00-00-01_00:00:00')
    parser.add_argument('-o', dest="outputInterval", default='00-00-01_00:00:00')
    args = parser.parse_args()

    run_testcase(args.nProcs,
                 args.runDuration,
                 args.outputInterval)
