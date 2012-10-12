#!/bin/env python
#############################################################################
# reskim.py - reskim files on dCache
#
# (C) RWTH Aachen University III. Physikalisches Insitut A 2012
# Author: M. Weber

import os
import sys
import optparse
import ConfigParser
import subprocess
import shutil
import getpass
import time

######################################################################
# Configuration

# The storage element (protocol://host:port)
se="srm://grid-srm.physik.rwth-aachen.de:8443"

# The dcache server (protocol://host)
dcap="dcap://grid-dcap.physik.rwth-aachen.de"

# The grid ftp server (only the host)
ftp="grid-ftp.physik.rwth-aachen.de"

# The path to your user area on the grid storage element (note the leading slash)
sp="/pnfs/physik.rwth-aachen.de/cms/store/user/"

# CONDOR template file for job submission
condor_template='''# Condor universe
universe     = vanilla

# Condor job description template
executable   = %(executable)s
arguments    = %(arguments)s
transfer_executable = True

transfer_input_files = %(inputfiles)s
transfer_output_files = %(outputfiles)s
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

# logging 
error   = %(stderr)s
output  = %(stdout)s
log     = %(log)s

# only send notifications on error
notification = Error

# copy environment variables from submitting machine, neccessary for ROOT
getenv = True

# which machines to prefer, rank on MIPS
rank = TARGET.Mips
requirements = (Mips > 5000 && KeyboardIdle>10000)

# queue this job
Queue
'''

# job starter template
job_starter='''#!/bin/bash
HOST=`uname -n`
echo "job_starter running on host $HOST"
echo "Current working dir is $PWD"
echo "Available disk space:"
df -h .
echo "Current directory contents:"
ls -l
#(( SLEEP= $RANDOM % 15 ))
#echo "Sleeping for $SLEEP seconds..."
#sleep $SLEEP
#echo "LD_LIBRARY_PATH contents:"
#echo $LD_LIBRARY_PATH | tr ':' '\\n'
#echo "Environment:"
#env
EXE=$1
shift
ARGS=$*
echo "Executable is $EXE"
echo "Arguments are $ARGS"
echo "Starting executable now..."
$EXE $ARGS
RC=$?
echo "Executable finished with return code $RC"
echo "Current directory contents:"
ls -l
echo "job_starter done."
exit $RC
'''

######################################################################
# functions

def getCommandOutput2(command):
    """Execute a command and get its output"""
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise RuntimeError, '%s failed with exit code %d' % (command, err)
    return data

def uberftpls(path, verbose=False):
    """List files on dCache"""
    fullpath = ftp+path
    files = []
    if verbose:
        print 'Listing ', fullpath
    # create temporary file with commands for uberftp
    commands = """cd %s
ls
""" % ( path )
    # call program and record both stderr and stdout
    p = subprocess.Popen(["uberftp", ftp], 
                         stdin=subprocess.PIPE, 
                         stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(commands)
    if len(stderr) > 0:
        print "Errors found:", stderr
        sys.exit(1)
    for line in stdout.splitlines():
        # remove uberftp prompt
        line = line.replace("uberftp> ", "")
        # split into file size and file name
        try:
            size = int(line.split()[4])
            fname = line.split()[8]
        except ValueError:
            continue
        except IndexError:
            continue
        if verbose:
            print size, fname
        files.append([size, fname])
    return files

def grep_file(filename, tokens, verbose=True):
    found = False
    searchfor = "|".join(tokens)
    child = os.popen('/bin/bash -c \'grep -Ei \"' + searchfor + '\" ' + filename + '\'')
    data = child.read().splitlines()
    rc = child.close()
    if len(data) > 0:
        found = True
        if verbose:
            print "Tokens matching %s found in file %s: " % ( searchfor, filename )
            for line in data:
                print line
    if rc == 2:
        raise Exception("Error executing grep \"query\" % for file %s " % (searchfor, filenam))
    return found

def check_log(logfilename, errors, warnings, requirements):
    """Check the log files and return True on success.
       errors is an iterable object containing strings marking errors in log file
       warnings is an iterable object containing strings marking warnings in log file
       requirements is an iterable object containing strings the log file that must appear

       If errors are found, print errors and return False
       If warnings are found, print warnings and return True
       If requirements are not found, print message and return False
    """

    # check for errors
    if grep_file(logfilename, errors):
        return False

    # check for warnings
    grep_file(logfilename, warnings)

    # check for requirements
    if requirements != None:
        if not grep_file(logfilename, requirements, verbose=False):
            print "Requirements not found in logfile %s: " % ( logfilename )
            return False

    # inform user
    return True

def setupInbox():
    """set up local files for transfer with CONDOR"""
    global options
    global inbox
    inbox = []

    # binary files
    myDir = options.outputdir
    os.system("mkdir -p " + myDir)

    # copy executable in place
    src = 'reskim'
    dest = myDir + '/reskim'
    shutil.copy(src, dest)
    inbox.append(dest)

    # create job starter file
    job_starter_filename = options.outputdir+"/job_starter.sh"
    exefile = open(job_starter_filename, "w")
    exefile.write(job_starter)
    exefile.close()
    os.system("chmod 755 " + job_starter_filename)

def wait_for_jobs(njobs):
    wait = True
    while wait:
        # check how many jobs are running already
        ntrial = 0
        while (ntrial < 3):
            try:
                joblist = getCommandOutput2("condor_q");
                break
            except RuntimeError as message:
                print message
                print "Waiting for condor to settle..."
                time.sleep(60)
                ntrial += 1
        jobcount = 0
        for line in joblist.splitlines():
            if getpass.getuser() in line:
                jobcount += 1
        # if the maximal job count is exceeded, sleep for a while
        if jobcount >= njobs:
            print "Reached " + str(njobs) + " jobs, sleeping for 30 s..."
            time.sleep(30)
        else:
            wait = False

def submit_condor_job(executable, arguments, inbox, outbox, jobname):
    global options
    stderr = "_".join((jobname, "stderr.log"))
    stdout = "_".join((jobname, "stdout.log"))
    cfgFile = "_".join((jobname, "condor.cfg"))
    repMap = {}
    repMap["executable"] = 'job_starter.sh'
    repMap["arguments"] = executable + ' ' + arguments
    repMap["inputfiles"] = " ".join(inbox)
    repMap["stdout"] = stdout
    repMap["stderr"] = stderr
    repMap["log"] = "condor.log"
    repMap["outputfiles"] = " ".join(outbox)
    content = condor_template % repMap
    jobfile = open(cfgFile, "w")
    jobfile.write(content)
    jobfile.close()
    wait_for_jobs(options.njobs)
    rc = getCommandOutput2("condor_submit " + cfgFile)
    return rc

######################################################################
# main
def main():
    usage = """usage: %prog [options] dcache-dir output-dir

Will create a and submit jobs from the given dcache directory to CONDOR.  The
output files will be stored in the output-dir, please use a local disk and not
a network disk."""
    optParser = optparse.OptionParser(usage)
    defaulttemplate="condor_template.cfg"
    defaultnjobs=20
    defaultmerge=20
    defaultanalysiscfg='default'
    optParser.add_option("-n", "--njobs", dest="njobs",
                         help="how many CONDOR jobs to run at once",
                         default=defaultnjobs)
    optParser.add_option("-m", "--merge", dest="nmerge",
                         help="split processing into this number of jobs",
                         default=defaultmerge)

    global options
    (options, args) = optParser.parse_args()
    if len(args) != 2:
        optParser.print_help()
        return 1

    options.inputdir = sp + '/' + getpass.getuser() + '/' +  args[0]
    options.outputdir = args[1]

    if os.environ['CMSSW_BASE'] == '':
        raise "You must setup correct CMSSW version for reskim to work"

    # get number of splits
    try:
        options.nmerge = int(options.nmerge)
    except ValueError:
        print "Merge has to be a number"
        return 1

    if options.nmerge <= 0:
        print "Merge must be greater than zero"
        return 2

    # get number of jobs
    try:
        options.njobs = int(options.njobs)
    except ValueError:
        print "Number of jobs has to be a number"
        return 1

    if options.njobs <= 0:
        print "Number of jobs must be greater than zero"
        return 2

    # setup local directory with all files that need to be shipped to the batch host
    setupInbox()

    print 'Adding files from directory', dcap+options.inputdir
    files = uberftpls(options.inputdir)
    filelist = []
    for (size, filename) in files:
        filelist.append(dcap + options.inputdir + '/' + filename)

    # need to go into this directory - submit job from there, output will get here
    os.system("mkdir -p " + options.outputdir)
    os.chdir(options.outputdir)

    # determine file splitting
    nfiles = len(filelist)
    print 'Found', nfiles, 'files'
    counter = 0
    while counter*options.nmerge < nfiles:
        start = counter*options.nmerge
        end = min((counter+1)*options.nmerge, nfiles)

        print "Job #", counter, ": ", end-start, " files, ", start, " - " , end-1

        executable = './reskim'
        jobName = 'reskim_' + str(counter) 
        outputFile = jobName + '.root'
        outbox = [ outputFile ]
        arguments = outputFile + ' ' + " ".join(filelist[start:end])

        submit_condor_job(executable, arguments, 
                          inbox, outbox, jobName)

        counter = counter + 1

######################################################################
# Call main if not included in other script

if __name__ == "__main__":
    sys.exit(main())
