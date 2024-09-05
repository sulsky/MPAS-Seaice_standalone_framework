import xml.etree.ElementTree as ET



columnPools = [
    "tracers",
    "tracers_aggregate",
    "icestate",
    "atmos_coupling",
    "atmos_forcing",
    "alternative_atmos_forcing",
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
                      "description"]


filenamein = "Registry.xml"

tree = ET.parse(filenamein)
root = tree.getroot()

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
            print("  ", var.tag, var.attrib["name"])
            print("     ", var.attrib["dimensions"])
            print("     ", var.attrib["dimensions"].replace("nCells","nParticles"))

            var_out = ET.SubElement(var_struct_out, "var")
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

outString = ET.tostring(out,encoding="unicode")

for attribute in attributesToFormat:
    outString = outString.replace(" "+attribute+"=","\n         "+attribute+"=")

outString = outString.replace(" />","\n    />")
outString = outString.replace("</var_struct>","</var_struct>\n")
outString = outString.replace("packuges","packages")

print(outString)

fileout = open("Registry_mpm_column.xml","w")
fileout.write(outString)
fileout.close()
