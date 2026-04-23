import os

#-------------------------------------------------------------------------------

def setup_testcase():

    MPAS_SEAICE_DOMAINS_DIR = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
    if (MPAS_SEAICE_DOMAINS_DIR is None):
        raise Exception("MPAS_SEAICE_DOMAINS_DIR must be set")

    domainScript = MPAS_SEAICE_DOMAINS_DIR + "/domain_QU30km_icebergs_antarctica/get_domain.py"

    os.system("python %s" %(domainScript))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    setup_testcase()
