import argparse
from netCDF4 import Dataset
import os

#--------------------------------------------------------------------

def add_deldyn_to_ics(dynamicsTimeStep):

    reses = ["2562","10242","40962","163842"]

    icTypes = ["cosine_bell","slotted_cylinder"]

    for icType in icTypes:

        print("icType: ", icType)

        for res in reses:

            print("  Res: ", res)

            icFilename = "ic_%s_%s.nc" %(icType, res)

            if (not os.path.isfile(icFilename)):
                raise Exception("IC file missing: "+icFilename)

            filein = Dataset(icFilename,"a")

            uVelocity = filein.variables["uVelocity"][:]
            vVelocity = filein.variables["vVelocity"][:]

            try:
                filein.createDimension("TWO",2)
            except:
                pass

            try:
                deluDyn = filein.createVariable("deluDyn","d",dimensions=["nVertices","TWO"])
            except:
                deluDyn = filein.variables["deluDyn"][:]
            deluDyn[:,0] = uVelocity[:] * dynamicsTimeStep
            deluDyn[:,1] = vVelocity[:] * dynamicsTimeStep

            filein.close()

#--------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', dest='dynamicsTimeStep', type=float)

    args = parser.parse_args()

    add_deldyn_to_ics(args.dynamicsTimeStep)
