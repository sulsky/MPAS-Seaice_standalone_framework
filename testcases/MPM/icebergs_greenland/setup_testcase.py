import os

#-------------------------------------------------------------------------------

def setup_testcase():

    MPAS_SEAICE_DOMAINS_DIR = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
    if (MPAS_SEAICE_DOMAINS_DIR is None):
        raise Exception("MPAS_SEAICE_DOMAINS_DIR must be set")

    domainScript = MPAS_SEAICE_DOMAINS_DIR + "/domain_QU30km_icebergs_greenland/get_domain.py"

    os.system("python %s" %(domainScript))

    os.symlink("../icebergs_antarctica/namelist.seaice.wagner_momentum.noseaice","namelist.seaice")
    os.symlink("../icebergs_antarctica/streams.seaice","streams.seaice")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    setup_testcase()
