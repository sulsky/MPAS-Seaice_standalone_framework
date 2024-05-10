import os, shutil
from compare_mpas_files import compare_files
from testing_utils import *

#-------------------------------------------------------------------------

def restartability(mpasDevelopmentDir,
                   SAFrameDirDev,
                   domainsDir,
                   domain,
                   configuration,
                   options,
                   check,
                   oversubscribe,
                   np1,
                   np2):

    # test time parameters
    #testTime  = "12:00:00"
    #testTime2 = "24:00:00"
    #restart_file = "restart.2000-01-02_00.00.00.nc"
    testTime  = "02:00:00"
    testTime2 = "04:00:00"
    restart_file = "restart.2000-01-01_04.00.00.nc"

    # find available directory name
    iTest = 1
    dirExists = True
    while (dirExists):
        testDir = "restartability_%i.%s.%s" %(iTest,configuration,domain)
        iTest = iTest + 1
        dirExists = os.path.isdir(testDir)

    # make a test directory
    create_test_directory(testDir)

    title = "Test: Restartability, Configuration: %s, Domain: %s" %(configuration,domain)

    print_colour(title, "title")

    logfile = open("log_test.txt","w")
    logfile.write(title)

    # base run
    nProcs = np1

    nmlChanges = {"seaice_model": {"config_run_duration":testTime2}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":testTime}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("base",
                  mpasDevelopmentDir,
                  SAFrameDirDev,
                  domainsDir,
                  domain,
                  configuration,
                  nmlChanges,
                  streamChanges,
                  nProcs,
                  logfile,
                  oversubscribe) != 0):
        run_failed("restartability")
        os.chdir("..")
        return 1

    # first restart run
    nProcs = np1

    nmlChanges = {"seaice_model": {"config_run_duration":testTime}}
    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = [{"streamName":"restart", "attributeName":"output_interval", "newValue":testTime}, \
                     {"streamName":"output" , "attributeName":"output_interval", "newValue":"none"}]

    if (run_model("restart1",
                  mpasDevelopmentDir,
                  SAFrameDirDev,
                  domainsDir,
                  domain,
                  configuration,
                  nmlChanges,
                  streamChanges,
                  nProcs,
                  logfile,
                  oversubscribe) != 0):
        run_failed("restartability")
        os.chdir("..")
        return 1

    # restart
    nProcs = np2

    bgcRestart = False
    if ("bgc" in options.keys() and options["bgc"] == "True"):
        bgcRestart = True

    snowModsRestart = False
    if ("snow_tracer_physics" in options.keys() and options["snow_tracer_physics"] == "True"):
        snowModsRestart = True

    if (not bgcRestart):
        if (not snowModsRestart):
             nmlChanges = {"seaice_model": {"config_start_time":"file"},
                           "restart": {"config_do_restart":True}}
        else:
             nmlChanges = {"seaice_model": {"config_start_time":"file"},
                      "restart": {"config_do_restart":True,
                                  "config_do_restart_snow_density":True,
                                  "config_do_restart_snow_grain_radius":True}}
    else:
        nmlChanges = {"seaice_model": {"config_start_time":"file"},
                      "restart": {"config_do_restart":True,
                                  "config_do_restart_bgc":True,
                                  "config_do_restart_hbrine":True}}

    if (check):
        nmlChanges["unit_test"] = {"config_testing_system_test":True}
    nmlChanges = add_pio_namelist_changes(nmlChanges, nProcs)

    streamChanges = []

    os.system("cp -rf restart1 restart2")
    if (restart_model("restart2",
                      nmlChanges,
                      streamChanges,
                      nProcs,
                      logfile,
                      oversubscribe,
                      ["particles.nc"]) != 0):
        run_failed("restartability")
        os.chdir("..")
        return 1


    # compare
    file1 = "./base/restarts/%s" %(restart_file)
    file2 = "./restart2/restarts/%s" %(restart_file)
    logfile.write("file1: %s\n" %(file1))
    logfile.write("file2: %s\n" %(file2))

    ignoreVarname = ["cellsOnCell","verticesOnCell","edgesOnEdge","edgesOnCell"]
    if (check):
        ignoreVarname.append("testArrayParallelism")
        ignoreVarname.append("testArrayReproducibility")

    nErrorsArray, nErrorsNonArray = compare_files(file1,file2,logfile,ignoreVarname)

    failed = test_summary(nErrorsNonArray, nErrorsArray, logfile, "restartability")

    logfile.close()

    os.chdir("..")

    return failed

#-------------------------------------------------------------------------
