#!/bin/bash
######################################################################
# This script manages the preparation of CMSSW skimming jobs         #
#                                                                    #
# (C) 2011 Martin Weber                                              #
#                                                                    #
######################################################################
SKIM_VERSION="v7-6"

function usage()
{
    NAME=`basename $0`
cat <<EOF
Preparation of CMSSW skimming jobs, either MC or data

SYNOPSIS: $NAME [ -d | -m ] 
                -u username -s dataset -v version -t tag -g globaltag 
                [ -j json-file -n njobs | -low pthatlow -high pthathigh ]

OPTIONS:
  -h               This help message
  -d               Run on data
  -m               Run on Monte-Carlo
  -u username      Your grid user name, e.g. mweber
  -s dataset       Dataset path, e.g. /DoubleMu/Run2011A-May10ReReco-v1/AOD
  -v version       Any version tag you like, defaults to skimmer version ($VERSION)
  -t tag           Any tag that allows you to identify this sample, e.g. "data"
  -g globaltag     The global tag to be used, e.g. "GR_R_42_V21A::All"
  -f               Force overwriting of existing directories or files

Only for MC:
  -l pthatlow      A value for specifying a lower cut on pthat
  -i pthathigh     A value for specifying an upper cut on pthat

Only for DATA:
  -j json-file     The json file to be used for good qualified runs
  -n njobs         Number of jobs to submit (default: $NJOBS)
EOF
}

function check_file()
{
    if [[ ! -f $1 ]] ; then
	echo $2
	exit 1
    fi

}

function createconfigfiles()
{
    # create CMSSW configuration file from template
    echo "Creating CMSSW configuration file..."
    if [[ -n $DATA ]] ; then
	ISDATA=True
	DATA=
	MC="#"
        # copy JSON to local directory
	cp $JSON $JOBDIR/
    else
	ISDATA=False
	DATA="#"
	MC=
    fi
    PYFILE="$JOBDIR/cmssw.py"
    sed -e s/@ISDATA@/$ISDATA/g -e s/@GLOBALTAG@/$GLOBALTAG/g $CMSSWCFG >> $PYFILE

    # create CRAB configuration file from template
    echo "Creating crab configuration file..."
    sed -e s:@DATASET@:$DATASET:g -e s:@DATA@:$DATA:g -e s:@MC@:$MC:g \
	-e s:@EMAIL@:$USER@cern.ch:g -e s:@REMOTEDIR@:$REMOTE_DIR:g \
        -e s:@JSON@:$JSON:g -e s:@NJOBS@:$NJOBS:g \
	$CRABCFG >> $JOBDIR/crab.cfg
}

function main()
{
    # initialize variables
    CMSSWCFG="cmssw_template.py"
    CRABCFG="crab_template.cfg"
    VERSION=${SKIM_VERSION}
    DATA=
    MC=
    USER=
    DATASET=
    TAG=
    GLOBALTAG=
    PTHAHTLOW=
    PTHATHIGH=
    JSON=
    NJOBS=500
    FORCE=
    # process command line options
    while getopts "hdmu:s:v:t:g:l:i:j:n:f" OPTION
    do
         case $OPTION in
             h)
                 usage
                 exit
                 ;;
             d)
		 DATA=1
		 ;;
	     m)
		 MC=1
		 ;;
	     u)
		 USER=$OPTARG
		 ;;
	     s)
		 DATASET=$OPTARG
		 ;;
	     v)
		 VERSION=$OPTARG
		 ;;
	     t)
		 TAG=$OPTARG
		 ;;
	     g)
		 GLOBALTAG=$OPTARG
		 ;;
	     l)
		 PTHATLOW=$OPTARG
		 ;;
	     i)
		 PTHATHIGH=$OPTARG
		 ;;
	     j)
		 JSON=$OPTARG
		 ;;
	     n)
		 NJOBS=$OPTARG
		 ;;
	     f)  FORCE=1
		 ;;
             ?)
                 usage
                 exit
                 ;;
         esac
    done
    # if no parameters given, give help
    if [[ $# == 0 ]] ; then
	usage
	exit 0
    fi
    # check parameter logic
    EXIT=
    if [[ -z $DATA && -z $MC ]] ; then
	echo "You must specify either -d or -m switch"
	EXIT=1
    fi
    if [[ -z $DATA && -n $NJOBS ]] || [[ -z $DATA && -n $JSON ]] ; then
	echo "-j and -n can only be used on data"
	EXIT=1
    fi
    if [[ -z $MC && -n $PTHATHIGH ]] || [[ -z $DATA && -n $PTHATLOW ]] ; then
	echo "-l and -i can only be used on MC"
	EXIT=1
    fi
    if [[ -z $USER || -z $DATASET || -z $VERSION || -z $TAG || -z $GLOBALTAG ]] ; then
	echo "You must specify the -u, -s, -t, and -g switches"
	EXIT=1
    fi
    if [[ -n $DATA && -z $JSON ]] ; then
	echo "When running on data, you must specify the JSON file with -j"
	EXIT=1
    fi
    if [[ -n $EXIT ]] ; then
	echo "Use -h to get help"
	exit 1
    fi
    if [[ -n $PTHATLOW || -n $PTHATHIGH ]] ; then
	echo "pthatlow and pthathight not yet implemented!"
	exit 2
    fi

    # check files
    check_file $CMSSWCFG "ERROR: Python configuration template $CMSSWCFG does not exist"
    check_file $CRABCFG "ERROR: Crab configuration template $CRABCFG does not exist"

    check_file $JSON "ERROR: JSON file $JSON does not exist"

    # Create working directory
    JOBDIR="CRAB-$TAG-$VERSION"
    if [[ -d $JOBDIR ]]; then
	if [[ -z $FORCE ]] ; then 
	    echo "ERROR: Working directory $JOBDIR already exists"
	    exit 4
	else
	    echo "Deleting existing job directory $JOBDIR"
	    rm -rf $JOBDIR
	fi
    fi
    echo "Creating job directory $JOBDIR"
    mkdir $JOBDIR

    # Create output directory
    STORAGE_SERVER="srm://grid-srm.physik.rwth-aachen.de:8443"
    REMOTE_DIR="output/$VERSION/$TAG"
    STORAGE_DIR="pnfs/physik.rwth-aachen.de/cms/store/user/$USER/$REMOTE_DIR"
    OUTPUT_DIR="$STORAGE_SERVER/$STORAGE_DIR"
    srmls $OUTPUT_DIR >& /dev/null
    if [[ $? == 0 ]] ; then
	if [[ -z $FORCE ]]; then
	    echo "Directory $REMOTE_DIR does exist, exiting..."
	    exit 2
	else
	    echo "Deleting existing ouput directory $REMOTE_DIR"
	    srmrmdir -recursive $OUTPUT_DIR
	fi
    fi
    echo "Creating output directory $STORAGE_DIR"
    srmmkdir $OUTPUT_DIR >& /dev/null

    # create configuration files
    createconfigfiles

    cat <<EOF
Done. You may proceed with copy & paste of these commands:

voms-proxy-init -voms cms:/cms/dcms -valid 164:00
cd $JOBDIR
crab -create
crab -submit
crab -status
EOF
}

# Call main 
main $@
