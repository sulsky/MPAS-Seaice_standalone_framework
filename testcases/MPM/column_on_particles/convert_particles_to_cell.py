from netCDF4 import Dataset
import argparse

#-------------------------------------------------------------------------------

# converts MPAS-Seaice-MPM particle fields to be relabeled as cell fields

def convert_particles_to_cell(filenameIn,
                              filenameOut):

    filein = Dataset(filenameIn,"r")
    fileout = Dataset(filenameOut,"w")

    nParticles = len(filein.dimensions["nParticles"])

    newDimensions = []

    # dimensions
    for variable in filein.variables:

        var = filein.variables[variable]
        if ("nParticles" in var.dimensions):
            dimensions = var.dimensions
            for dimension in dimensions:
                if (dimension not in newDimensions):
                    if (dimension == "nParticles"):
                        fileout.createDimension("nCells",nParticles)
                        newDimensions.append("nParticles")
                    else:
                        fileout.createDimension(dimension,len(filein.dimensions[dimension]))
                        newDimensions.append(dimension)

        # include xtime
        if (variable == "xtime"):
            dimensions = var.dimensions
            for dimension in dimensions:
                if (dimension not in newDimensions):
                    fileout.createDimension(dimension,len(filein.dimensions[dimension]))
                    newDimensions.append(dimension)

    # variables
    for variable in filein.variables:

        # dimensions in variable
        var = filein.variables[variable]
        dimensions = var.dimensions
        newDimensions = []
        for dimension in dimensions:
            if (dimension == "nParticles"):
                newDimensions.append("nCells")
            else:
                newDimensions.append(dimension)

        if ("nParticles" in var.dimensions):
            newVar = fileout.createVariable(variable.rstrip("MP"),"d",newDimensions)
            newVar[:] = var[:]
            for attrName in var.ncattrs():
                setattr(newVar, attrName, getattr(var, attrName))

        if (variable == "xtime"):
            newVar = fileout.createVariable(variable,var.datatype, var.dimensions)
            newVar[:] = var[:]
            for attrName in var.ncattrs():
                setattr(newVar, attrName, getattr(var, attrName))

    filein.close()
    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest="filenameIn", required=True)
    parser.add_argument('-o', dest="filenameOut", required=True)
    args = parser.parse_args()

    convert_particles_to_cell(args.filenameIn,
                              args.filenameOut)
