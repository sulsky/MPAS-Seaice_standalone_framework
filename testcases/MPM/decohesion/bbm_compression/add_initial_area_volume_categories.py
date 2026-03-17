import argparse
import numpy as np
import netCDF4 as nc
import os

#--------------------------------------------------------------------

def add_initial_area_volume_categories():

   icFilename = 'particles.nc'
   if (not os.path.isfile(icFilename)):
      raise Exception("Mising particle file: "+icFilenam)

   cmd = "mv particles.nc particle_sourcefile.nc"
   os.system(cmd)

   # make a new netcdf file with proper category variables
   ds_orig = nc.Dataset('particle_sourcefile.nc', 'r')
   ds_new = nc.Dataset(icFilename, 'w', format='NETCDF3_CLASSIC')

   # create new dimension
   nCategories = 5
   ds_new.createDimension('nCategories', nCategories)
   # copy other dimensions
   for name, dimension in ds_orig.dimensions.items():
      if name != 'nCategories':
          ds_new.createDimension(name, dimension.size if not dimension.isunlimited() else None)

   # copy other (non-category) variables from the old file to the new file
   exclude_from_copy = ["iceAreaCategoryMP", "iceVolumeCategoryMP", "iceAreaCellMP", "iceVolumeCellMP"]
   for name, variable in ds_orig.variables.items():
      if (name not in exclude_from_copy):
         var = ds_new.createVariable(name, variable.dtype, variable.dimensions)
         var[:] = variable[:]
         for attr_name, attr_value in variable.__dict__.items():
            if attr_name not in ['_FillValue']:
               setattr(var, attr_name, attr_value)

   # copy attributes
   for attr_name, attr_value in ds_orig.__dict__.items():
      setattr(ds_new, attr_name, attr_value)

   # create the category variables and data
   nParticles = len(ds_orig.dimensions["nParticles"])
   iceAreaCategoryMP =   np.zeros((nParticles,nCategories))
   iceVolumeCategoryMP = np.zeros((nParticles,nCategories))

   iceAreaCellMP =   np.zeros(nParticles)
   iceVolumeCellMP = np.zeros(nParticles)

   categoryThicknessLimits = \
               [0, 0.6, 1.4, 2.4, 3.6, 100000000]

   iceAreaCategoryInit = np.zeros(nCategories)
   iceAreaCategoryInit[0] = 0.05
   iceAreaCategoryInit[1] = 0.1
   iceAreaCategoryInit[2] = 0.3
   iceAreaCategoryInit[3] = 0.35
   iceAreaCategoryInit[4] = 0.2

   iceThicknesses = np.zeros(nCategories)
   for iCategory in range(0,nCategories-1):
       iceThicknesses[iCategory] = 0.5 * (categoryThicknessLimits[iCategory] +
                                          categoryThicknessLimits[iCategory+1])
       iceThicknesses[nCategories-1] = categoryThicknessLimits[nCategories-1] + 1.0

   for iParticle in range(0, nParticles):
       for iCategory in range(0, nCategories):
          iceAreaCategoryMP[iParticle,iCategory] = iceAreaCategoryInit[iCategory]
          iceAreaCellMP[iParticle] = np.sum(iceAreaCategoryMP[iParticle,:])

       for iCategory in range(0, nCategories):
           iceVolumeCategoryMP[iParticle,iCategory] = iceThicknesses[iCategory] \
                                                      * iceAreaCategoryMP[iParticle,iCategory]
           iceVolumeCellMP[iParticle] = np.sum(iceVolumeCategoryMP[iParticle,:])

   # add category data to new file
   var = ds_new.createVariable("iceAreaCategoryMP", "d", dimensions=["nParticles","nCategories","ONE"])
   var[:,:,0] = iceAreaCategoryMP[:,:]

   var = ds_new.createVariable("iceVolumeCategoryMP", "d", dimensions=["nParticles","nCategories","ONE"])
   var[:,:,0] = iceVolumeCategoryMP[:,:]

   var = ds_new.createVariable("iceAreaCellMP", "d", dimensions=["nParticles"])
   var[:] = iceAreaCellMP[:]

   var = ds_new.createVariable("iceVolumeCellMP", "d", dimensions=["nParticles"])
   var[:] = iceVolumeCellMP[:]

   ds_new.close()
   ds_orig.close()
   os.remove('particle_sourcefile.nc')

#--------------------------------------------------------------------

if __name__ == "__main__":

    add_initial_area_volume_categories()
