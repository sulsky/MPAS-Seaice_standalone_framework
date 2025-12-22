from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def create_iceberg_categories():

    nIcebergCategories = 10

    icebergCategoryLength = [60.0,100.0,200.0,350.0,500.0,700.0,900.0,1200.0,1600.0,2200.0]
    icebergCategoryHeight = [40.0,67.0,133.0,175.0,250.0,250.0,250.0,250.0,250.0,250.0]
    icebergCategoryFluxScaling = [2000.0,200.0,50.0,20.0,10.0,5.0,2.0,1.0,1.0,1.0]
    icebergCategoryCalvingFraction = [0.25,0.12,0.15,0.18,0.12,0.07,0.03,0.03,0.03,0.02]

    icebergCategoryLength = np.array(icebergCategoryLength)
    icebergCategoryHeight = np.array(icebergCategoryHeight)
    icebergCategoryFluxScaling = np.array(icebergCategoryFluxScaling)
    icebergCategoryCalvingFraction = np.array(icebergCategoryCalvingFraction)

    fileout = Dataset("iceberg_categories.nc","w",format="NETCDF3_CLASSIC")

    fileout.src = "T. Martin and A. Adcroft (2010), Parameterizing the fresh-water flux from land ice to ocean with interactive icebergs in a coupled climate model, Ocean Modelling, 34, 111-124, https://doi.org/10.1016/j.ocemod.2010.05.001"

    fileout.createDimension("nIcebergCategories",nIcebergCategories)

    var = fileout.createVariable("icebergCategoryLength","d", dimensions=["nIcebergCategories"])
    var[:] = icebergCategoryLength[:]

    var = fileout.createVariable("icebergCategoryHeight","d", dimensions=["nIcebergCategories"])
    var[:] = icebergCategoryHeight[:]

    var = fileout.createVariable("icebergCategoryFluxScaling","d", dimensions=["nIcebergCategories"])
    var[:] = icebergCategoryFluxScaling[:]

    var = fileout.createVariable("icebergCategoryCalvingFraction","d", dimensions=["nIcebergCategories"])
    var[:] = icebergCategoryCalvingFraction[:]

    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_iceberg_categories()
