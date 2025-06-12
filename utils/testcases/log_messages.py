import sys

#-------------------------------------------------------------------------------

def log_message(message,
                color,
                doPrint=True,
                logFile=None):

    try:
        import colorama
        if (color == "red"):
            coloramaColor = colorama.Fore.RED
        elif (color == "green"):
            coloramaColor = colorama.Fore.GREEN
        elif (color == "magenta"):
            coloramaColor = colorama.Fore.MAGENTA
        elif (color == "yellow"):
            coloramaColor = colorama.Fore.YELLOW
        else:
            raise Exception("Unsupported color")
        if (doPrint):
            print(coloramaColor + message + colorama.Style.RESET_ALL)
        if (logFile is not None):
            logFile.write(coloramaColor + message + colorama.Style.RESET_ALL + '\n')
    except ImportError:
        if (doPrint):
            print(message)
        if (logFile is not None):
            logFile.write(message1 + '\n')

    if (doPrint):
        sys.stdout.flush()
    if (logFile is not None):
        logFile.flush()

#-------------------------------------------------------------------------------

def regression_summary(nErrorsNonArray, nErrorsArray, logFile, testname):

    if (nErrorsNonArray == 0):
        log_message("No non-contents errors for %s" %(testname), "green", False, logFile)
    else:
        log_message("%i non-content errors for %s" %(nErrorsNonArray, testname), "red", False, logFile)

    if (nErrorsArray == 0):
        log_message("PASS: Test %s passed" %(testname), "green", False, logFile)
    else:
        log_message("FAIL!: Test %s failed" %(testname), "red", False, logFile)

#-------------------------------------------------------------------------------
