import os

#-------------------------------------------------------------------------------

def get_testcase_data_square():

    MPAS_SEAICE_STANDALONE_DATA = os.environ.get('MPAS_SEAICE_STANDALONE_DATA')
    if (MPAS_SEAICE_STANDALONE_DATA is None):
        raise Exception("MPAS_SEAICE_STANDALONE_DATA must be set")

    filenames = ["ic_hex.nc",
                 "ic_quad.nc",
                 "square_mesh_80x80_culled.nc",
                 "square_mesh_82x94_culled.nc"]

    for filenameDst in filenames:

        filenameSrc = "%s/testcases/square/%s" %(MPAS_SEAICE_STANDALONE_DATA,filanameDst)

        if (not os.path.isfile(filanameDst)):
            os.symlink(filenameSrc, filanameDst)

    args = ["wget", dirName+filename]

    process = subprocess.Popen(args, stdout=subprocess.PIPE)

    while process.poll() is None:
        line = process.stdout.readline()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    get_testcase_data_square()
