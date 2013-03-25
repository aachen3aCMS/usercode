#!/bin/env python
import os
import sys
import inspect
import ConfigParser
import stat
import optparse

collectpy_dir = os.path.split(os.path.realpath(inspect.getfile(inspect.currentframe())))[0]
sys.path.append(collectpy_dir)

import reskim

######################################################################
# functions

def resubmit(condor_jobfile):
    global options
    myDir = os.path.dirname(condor_jobfile)
    os.chdir(myDir)
    reskim.wait_for_jobs(options.njobs)
    print "Resubmitting", condor_jobfile
    rc = reskim.getCommandOutput2("condor_submit " + condor_jobfile)
    return rc

def check(fname):
    """Check log files for errors and warnings"""

    good = True

    errors = [ "exception", "segmentation", "error", "err:", "return code 127"]
    warnings = [ "warning", "wrn" ]
    requirements = [ "Executable finished with return code 0" ]
    
    os.chdir(collectpy_dir)
    print "checking:",fname.replace("_condor.cfg", "_stdout.log")
    success = reskim.check_log(fname.replace("_condor.cfg", "_stdout.log"),
                               errors, warnings, requirements)
    #pdf_success = reskim.check_log(fname.replace("_condor.cfg", "_stdout.log"),
                                   #no_errors, warnings, requirements)
    if not success:
        good = False
        resubmit(fname)
    #else:
        #success = reskim.check_log(fname.replace("_condor.cfg", "_stderr.log"),
                                   #errors, warnings, None)
        #if not success:
            #good = False
            #resubmit(fname)
    return good

######################################################################
# main
def main():
    global options

    usage = """usage: %prog [outputdir]

Will check output of CONDOR jobs and resubmit failed jobs"""
    optParser = optparse.OptionParser(usage)
    defaultnjobs=20
    optParser.add_option("-n", "--njobs", dest="njobs",
                         help="how many jobs to run at once if resubmitting",
                         default=defaultnjobs)

    (options, args) = optParser.parse_args()
    if len(args) != 1:
        optParser.print_help()
        return 1

    try:
        ROOTSYS=os.environ['ROOTSYS']
    except KeyError:
        print("You must setup correct ROOTSYS for AdvancedROOTAnalyzer to work")
        return 1

    # get number of jobs
    try:
        options.njobs = int(options.njobs)
    except ValueError:
        print "Number of jobs has to be a number"
        return 1

    if options.njobs <= 0:
        print "Number of jobs must be greater than zero"
        return 2

    # setup local arguments
    options.selection = args[0]

    # list files in directory
    print "Reading local file names..."
    filelist = os.listdir(args[0])

    # loop over files and check them
    errors = False
    for fname in filelist:
        if fname.startswith("reskim_") and fname.endswith("_condor.cfg"):
            # configuration file found, check log files
            print "Job file found:", fname
            result = check(args[0].rstrip('/') + '/' + fname) 
            if result:
                print "Job OK!"
            errors = result and errors

    if errors:
        print "*******************************************************************************"
        print "* There were errors, check script output                                      *" 
        print "*******************************************************************************"
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
