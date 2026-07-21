import xml.etree.ElementTree as ET
import argparse
import os

#-------------------------------------------------------------------------------

def create_mpm_registry_file(filenameIn,
                             filenameOut):

    columnPools = [
        "tracers",
        "tracers_aggregate",
        "icestate",
        "atmos_coupling",
        "atmos_forcing",
        "alternative_atmos_forcing",
        "floe_size_distribution",
        "wave_coupling",
        "ocean_coupling",
        "general_column",
        "ridging",
        "melt_growth_rates",
        "atmos_fluxes",
        "ocean_fluxes",
        "ocean_atmosphere",
        "shortwave",
        "drag",
        "ponds",
        "snow",
        "aerosols",
        "biogeochemistry",
        "diagnostics_biogeochemistry",
        "initial",
        "diagnostics"]

    otherPools = [
        "velocity_solver",
        "boundary"]

    attributesToFormat = ["name_in_code",
                          "type",
                          "dimensions",
                          "packages",
                          "default_value",
                          "description"]

    tree = ET.parse(filenameIn)
    root = tree.getroot()

    varnamesOut = open("varnames.txt","w")

    out = ET.Element("registry")

    for var_struct in root:
        if (var_struct.tag == "var_struct" and var_struct.attrib["name"] in columnPools):

            var_struct_out = ET.SubElement(out, "var_struct")
            var_struct_out.attrib["name"] = var_struct.attrib["name"] + "_MPM"
            if ("time_levs" in var_struct.attrib):
                var_struct_out.attrib["time_levs"] = var_struct.attrib["time_levs"]
            if ("packages" in var_struct.attrib):
                packagesOut = ""
                for package in var_struct.attrib["packages"].split(";"):
                    packagesOut += package + "MP;"
                var_struct_out.attrib["packuges"] = packagesOut.rstrip(";")
            else:
                var_struct_out.attrib["packuges"] = "pkgMPMColumn"

            print(var_struct.attrib["name"], var_struct.attrib["name"]+"_MPM")

            for var in var_struct:
                #print("  ", var.tag, var.attrib["name"])
                #print("     ", var.attrib["dimensions"])
                #print("     ", var.attrib["dimensions"].replace("nCells","nParticles"))

                var_out = ET.SubElement(var_struct_out, "var")
                varnamesOut.write("%s\n" %(var.attrib["name"]))
                var_out.attrib["name"] = var.attrib["name"] + "MP"
                var_out.attrib["name_in_code"] = var.attrib["name"]
                var_out.attrib["type"] = var.attrib["type"]
                var_out.attrib["dimensions"] = var.attrib["dimensions"].replace("nCells","nParticles")
                if ("packages" in var.attrib):
                    packagesOut = ""
                    for package in var.attrib["packages"].split(";"):
                        packagesOut += package + "MP;"
                    var_out.attrib["packages"] = packagesOut.rstrip(";")
                if ("description" in var.attrib):
                    var_out.attrib["description"] = var.attrib["description"]

    ET.indent(out, space="  ", level=0)

    varnamesOut.close()

    outString = ET.tostring(out,encoding="unicode")

    # add SSS 34 default to mpm column registry file
    sstIn  = "<var name=\"seaSurfaceSalinityMP\" name_in_code=\"seaSurfaceSalinity\" type=\"real\" dimensions=\"nParticles Time\" description=\"Sea surface salinity\" />"
    sstOut = "<var name=\"seaSurfaceSalinityMP\" name_in_code=\"seaSurfaceSalinity\" type=\"real\" dimensions=\"nParticles Time\" default_value=\"34.0\" description=\"Sea surface salinity\" />"
    if (sstIn not in outString):
        raise Exception("Could not find seaSurfaceSalinity string in outString")
    outString = outString.replace(sstIn,sstOut)

    for attribute in attributesToFormat:
        outString = outString.replace(" "+attribute+"=","\n         "+attribute+"=")

    outString = outString.replace(" />","\n    />")
    outString = outString.replace("</var_struct>","</var_struct>\n")
    outString = outString.replace("packuges","packages")

    outString = outString.replace("<registry>\n","")
    outString = outString.replace("</registry>","")

    filenameStart = os.path.dirname(os.path.abspath(__file__))+"/registry_mpm_column_start.txt"
    fileStart = open(filenameStart,"r")
    startString = fileStart.read()
    fileStart.close()

    fileout = open(filenameOut,"w")
    fileout.write(startString)
    fileout.write(outString)
    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='filenameIn', required=True)
    parser.add_argument('-o', dest='filenameOut', required=True)

    args = parser.parse_args()

    create_mpm_registry_file(args.filenameIn,
                             args.filenameOut)
