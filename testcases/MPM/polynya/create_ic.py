from netCDF4 import Dataset

#-------------------------------------------------------------------------------

def create_ic():

    # mesh data
    filein = Dataset("grid.nc","r")

    nVertices = len(filein.dimensions["nVertices"])

    filein.close()

    # ice file
    uVelocity = 0.5
    vVelocity = 0.0

    # deluDyn
    dynamicsTimeStep = 3600.0
    deluDynu = uVelocity * dynamicsTimeStep
    deluDynv = vVelocity * dynamicsTimeStep

    fileout = Dataset("ic.nc","w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nVertices",nVertices)
    fileout.createDimension("TWO",2)

    var = fileout.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = uVelocity

    var = fileout.createVariable("vVelocity","d",dimensions=["nVertices"])
    var[:] = vVelocity

    var = fileout.createVariable("deluDyn","d",dimensions=["nVertices","TWO"])
    var[:,0] = deluDynu
    var[:,1] = deluDynv

    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ic()
