from netCDF4 import Dataset
import argparse
import glob
import sys
sys.path.append("../../../testing")
from testing_utils import print_colour

#-------------------------------------------------------------------------------

def check_particle_positions():

    try:
        import colorama
        colorama.init()
    except ImportError:
        pass

    print()
    operatorMethods = ["mpmvar","mpas"]
    message = "Check particles positions"
    print_colour(message, "title")

    filenames1 = sorted(glob.glob("./output_mpmvar/particles*"))
    filenames2 = sorted(glob.glob("./output_mpas/particles*"))

    print("filename for mpmvar: %s" %(filenames1[-1]))
    print("filename for mpas: %s"  %(filenames2[-1]))

    file1 = Dataset(filenames1[-1], "r")
    file2 = Dataset(filenames2[-1], "r")

    nParticles1 = len(file1.dimensions["nParticles"])
    nParticles2 = len(file2.dimensions["nParticles"])

    posnMP1 = file1.variables["posnMP"][0,:]
    posnMP2 = file2.variables["posnMP"][0,:]

    statusMP1 = file1.variables["statusMP"][0,:]
    statusMP2 = file2.variables["statusMP"][0,:]

    cellIDCreationMP1 = file1.variables["cellIDCreationMP"][0,:]
    cellIDCreationMP2 = file2.variables["cellIDCreationMP"][0,:]

    creationIndexMP1 = file1.variables["creationIndexMP"][0,:]
    creationIndexMP2 = file2.variables["creationIndexMP"][0,:]

    file1.close()
    file2.close()


    # check number of active particles
    nActive1 = sum(statusMP1)
    nActive2 = sum(statusMP2)
    if (nActive1 != nActive2):
        message1 = "Difference in number of active material points detected" 
        message2 = "%s has %i acitve mps; %s has %i active MPs" %("mpmvar", nActive1, "mpas", nActive2)
        try:
            import colorama
            print(colorama.Fore.RED + message1 + colorama.Style.RESET_ALL)
            print(colorama.Fore.RED + message2 + colorama.Style.RESET_ALL)
        except ImportError:
            print(message1)
            print(message2)
        print()
        raise Exception("Difference in number of active material points")


    cellIDToIndex2 = {}
    for iParticle2 in range(0,nParticles2):
        if (statusMP2[iParticle2] == 1):
            cellIDToIndex2[(cellIDCreationMP2[iParticle2], creationIndexMP2[iParticle2])] = iParticle2


    error = False
    for iParticle1 in range(0,nParticles1):

        if (statusMP1[iParticle1] == 1):

            iParticle2 = cellIDToIndex2[(cellIDCreationMP1[iParticle1], creationIndexMP1[iParticle1])]

            dx = posnMP2[iParticle2,0] - posnMP1[iParticle1,0]
            dy = posnMP2[iParticle2,1] - posnMP1[iParticle1,1]
            dz = posnMP2[iParticle2,2] - posnMP1[iParticle1,2]

            if (dx != 0.0):
                print("dx different: ", dx, iParticle1, iParticle2, cellIDCreationMP1[iParticle1])
                error = True
            if (dy != 0.0):
                print("dy different: ", dy, iParticle1, iParticle2, cellIDCreationMP1[iParticle1])
                error = True
            if (dz != 0.0):
                print("dz different: ", dz, iParticle1, iParticle2, cellIDCreationMP1[iParticle1])
                error = True

    if (error):
        message = "Error: Differences in particle position detected for mpmvar and mpas"
        try:
            import colorama
            print(colorama.Fore.RED + message + colorama.Style.RESET_ALL)
        except ImportError:
            print(message)
        print()
        raise Exception(message)

    message = "No differences detected in particle position for mpmvar and mpas"
    try:
        import colorama
        print(colorama.Fore.GREEN + message + colorama.Style.RESET_ALL)
    except ImportError:
        print(message)
    print()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    check_particle_positions()
