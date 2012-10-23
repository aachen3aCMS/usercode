#!/bin/bash
######################################################################
# This script manages the preparation of CMSSW skimming jobs         #
#                                                                    #
# (C) 2011 Martin Weber                                              #
#                                                                    #
######################################################################

function usage()
{
    NAME=`basename $0`
cat <<EOF
Preparation of CMSSW skimming jobs, either MC or data

SYNOPSIS: $NAME [ -d | -m ]
                -s dataset -t tag -c cutconfig
                [ -g globaltag -j json-file -n njobs ]

OPTIONS:
  -h               This help message
  -m 		   Running on MC
  -d               Running on data
  -s dataset       Dataset path, e.g. /DoubleMu/Run2011A-May10ReReco-v1/AOD
  -t tag           Any tag that allows you to identify this sample, e.g. "data"
  -g globaltag     global tag, overriding value from configuration file
  -f               Force overwriting of existing directories or files
  -c 		   Configfile for skimming cuts (default cuts.cfg)

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
    sed -e s/@ISDATA@/$ISDATA/g -e s/@GLOBALTAG@/$GLOBALTAG/g -e s/@QSCALE_LOW@/$QSCALELOW/g -e s/@QSCALE_HIGH@/$QSCALEHIGH/g -e s/@ELEPT@/$ELEPT/g -e s/@ELEETA@/$ELEETA/g -e s/@NELE@/$NELE/g -e s/@NMUO@/$NMUO/g -e s/@MUOPTFIRST@/$MUOPTFIRST/g -e s/@MUOPTOTHER@/$MUOPTOTHER/g -e s/@MUOETA@/$MUOETA/g -e s/@PFELEPT@/$PFELEPT/g -e s/@PFELEETA@/$PFELEETA/g -e s/@NPFELE@/$NPFELE/g -e s/@COMMONSKIM@/$COMMONSKIM/g -e s/@PHOPT@/$PHOPT/g -e s/@PHOETA@/$PHOETA/g -e s/@NPHO@/$NPHO/g -e s/@TAUPT@/$TAUPT/g -e s/@TAUETA@/$TAUETA/g -e s/@NTAU@/$NTAU/g -e s/@PFJETPT@/$PFJETPT/g -e s/@PFJETETA@/$PFJETETA/g -e s/@NPFJET@/$NPFJET/g -e s/@ELEPT@/$ELEPT/g -e s/@ELEETA@/$ELEETA/g -e s/@NELE@/$NELE/g -e s/@MET0@/$MET0/g -e s/@MET1@/$MET1/g -e s/@MET2@/$MET2/g -e s/@PFHTC@/$PFHTC/g -e s/@HTC@/$HTC/g -e s/@TRIGGERCONTAINS@/$TRIGGERCONTAINS/g -e s/@DOTAU@/$DOTAU/g -e s/@DOPFELE@/$DOPFELE/g -e s/@MUOMINV@/$MUOMINV/g -e s/@MUODMINV@/$MUODMINV/g  -e s/@PYTHIA8@/$PYTHIA8/g -e s/@SHERPA@/$SHERPA/g -e s/@MATCHALL@/$MATCHALL/g -e s/@SUSYPAR@/$SUSYPAR/g -e s/@ISPYTHIASHOWERED@/$ISPYTHIASHOWERED/g $CMSSWCFG >> $PYFILE
	#~ sed -e s/@ISDATA@/$ISDATA/g -e s/@GLOBALTAG@/$GLOBALTAG/g $CMSSWCFG >> $PYFILE
    # create CRAB configuration file from template
    echo "Creating crab configuration file..."
    sed -e s:@DATASET@:$DATASET:g -e s:@DATA@:$DATA:g -e s:@MC@:$MC:g \
	-e s:@EMAIL@:$USER@cern.ch:g -e s:@REMOTEDIR@:$REMOTE_DIR:g \
        -e s:@JSON@:$JSON:g -e s:@NJOBS@:$NJOBS:g \
	$CRABCFG >> $JOBDIR/crab.cfg
}
COUNTER=0
function srm_mkdir()
{
    let COUNTER++
    if [[ $COUNTER > 5 ]]; then
	echo "Maximum recursion level reached - stopping"
	echo "This probably means something is wrong with the storage path."
	exit 3
    fi
    srmmkdir $1 >& /dev/null
    RC=$?
    if [[ $RC != 0 ]] ; then
	# start recursion
	srm_mkdir `dirname $1`
	srmmkdir $1 >& /dev/null
	RC=$?
	if [[ $RC != 0 ]] ; then
	    echo "Making directories recursively failed - this should not happen."
	    exit 4
	fi
    fi
    echo "Creating $1"
}

function main()
{
    # initialize variables
    CMSSWCFG="cmssw_template_crab.py"
    CRABCFG="crab_template.cfg"
    CUTCFG="SkimmingCuts.cfg"
    VERSION=
    DATA=
    MC=
    USER=
    DATASET=
    TAG=
    GLOBALTAG=
    QSCALELOW=-1
    QSCALEHIGH=-1
    JSON=
    NJOBS=500
    FORCE=
    # process command line options
    while getopts "hdmu:s:v:t:g:l:c:i:j:n:f" OPTION
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
	     s)
		 DATASET=$OPTARG
		 ;;
	     t)
		 TAG=$OPTARG
		 ;;
	     g)
		 GLOBALTAG=$OPTARG
		 ;;
	     l)
		 QSCALELOW=$OPTARG
		 ;;
	     i)
		 QSCALEHIGH=$OPTARG
		 ;;
	     j)
		 JSON=$OPTARG
		 ;;
	     n)
		 NJOBS=$OPTARG
		 ;;
	     c)  
		 CUTCFG=$OPTARG
		 ;;
	     f)  
		 FORCE=1
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

    if [[ -z $DATA && -n $JSON ]] ; then
	echo "-j can only be used on data"
	EXIT=1
    fi
    #~ if [[ -z $MC && -n $QSCALEHIGH ]] || [[ -z $DATA && -n $QSCALELOW ]] ; then
	#~ echo "-l and -i can only be used on MC"
	#~ EXIT=1
    #~ fi

    if [[ -z $DATASET || -z $TAG ]] ; then
	echo "You must specify the -s and -t switches"
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
    #~ if [[ -n $PTHATLOW || -n $PTHATHIGH ]] ; then
	#~ echo "pthatlow and pthathigh not yet implemented!"
	#~ exit 2
    #~ fi

    echo "Sufficient information given on command line"

    # check files
    check_file $CMSSWCFG "ERROR: Python configuration template $CMSSWCFG does not exist"
    check_file $CRABCFG "ERROR: Crab configuration template $CRABCFG does not exist"
    check_file $CUTCFG "ERROR: Skimming configuration file $CUTCFG does not exist"

    echo "All files found"

    NTEST="0."
    NTEST2="0"
    NTEST4="1"

    # read skimmung cuts and configuration from config file
    source $CUTCFG

    echo "Cuts read from file"

    # sanity checks on skimming cuts. Gives warnings, does not prevent you from doing something you really want to.

    # sanity checks on skimming cuts
    if  [[ $NELE == $NTEST2  && $NMUO == $NTEST2  && $NPFELE == $NTEST2 && $NPHO == $NTEST2 && $NTAU == $NTEST2 && $NPFJET == $NTEST2 && $MET0 == $NTEST && $MET1 == $NTEST && $MET2 == $NTEST && $COMMONSKIM == 'False' ]] ; then
	echo "WARNING! You have specified no object to cut on!"
    fi
    if [[ $NELE == $NTEST2 &&  $ELEPT != $NTEST ]]; then
	echo "WARNING! elept > 0 but nele = 0. Cut has no effect!"
    fi
    if [[ $NMUO == $NTEST2  &&  $MUOPTFIRST != $NTEST ]]; then
	echo "WARNING! muopt > 0 but nmuo == 0. Cut has no effect!"
    fi
    if [[ $NPFELE == $NTEST2  &&  $PFELEPT != $NTEST ]]; then
	echo "WARNING! pfelept > 0 but npfele == 0. Cut has no effect!"
    fi
    if [[ $NTAU == $NTEST2  &&  $TAUPT != $NTEST ]]; then
	echo "WARNING! taupt > 0 but ntau == 0. Cut has no effect!"
    fi
    if [[ $NPHO == $NTEST2  &&  $PHOPT != $NTEST ]]; then
	echo "WARNING! phopt > 0 but npho == 0. Cut has no effect!"
    fi
    if [[ $NPFJET == $NTEST2  &&  $PFJETPT != $NTEST ]]; then
	echo "WARNING! pfjetpt > 0 but npfjet == 0. Cut has no effect!"
    fi
    if [[ $NELE != $NTEST2 ]]; then
	echo "Preparing Skim requiring $NELE electron(s) with pt > $ELEPT GeV and eta < $ELEETA"
    fi
    if [[ $NMUO != $NTEST2 ]]; then
	echo "Preparing Skim requiring $NMUO muon(s) with pt > $MUOPTFIRST GeV for the first and $MUOPTOTHER GeV for the others and eta < $MUOETA"
    fi
    if [[ $NPFELE != $NTEST2 ]]; then
	echo "Preparing Skim requiring $NPFELE PFelectron(s) with pt > $PFELEPT GeV and eta < $PFELEETA"
    fi
    if [[ $NPHO != $NTEST2 ]]; then
	echo "Preparing Skim requiring $NPHO photon(s) with pt > $PHOPT GeV and eta < $PHOETA"
    fi
    if [[ $NTAU != $NTEST2 ]]; then
	echo "Preparing Skim requiring $NTAU tau(s) with pt > $TAUPT GeV and eta < $TAUETA"
    fi
    if [[ $NPFJET != $NTEST2 ]]; then
	echo "Preparing Skim requiring $NPFJET pfjet(s) with pt > $PFJETPT GeV and eta < $PFJETETA"
    fi
    if [[ $MET0 != $NTEST ]]; then
	echo "Preparing Skim requiring MET0 > $MET0 GeV"
    fi
    if [[ $MET1 != $NTEST ]]; then
	echo "Preparing Skim requiring MET1 > $MET2 GeV"
    fi
    if [[ $MET2 != $NTEST ]]; then
	echo "Preparing Skim requiring MET2 > $MET2GeV"
    fi
    if [[ $COMMONSKIM != 'False' && $COMMONSKIM != 'True' ]] ; then
	echo "You have specified neither 'True' nor 'False' for COMMONSKIM, exiting"
	return 1
    fi
    if [[ $COMMONSKIM == 'True' ]] ; then 
	echo "======================================================================"
	echo "= Preparing common skim =============================================="
	echo "======================================================================"
    fi

    echo "Preparing crab task with following options:"
    NTEST3="False"

    if [[ $DOTAU == $NTEST3 ]]; then
	echo "Taus are OFF!"
    fi
    if [[ $DOTAU != $NTEST3 ]]; then
	echo "Taus are ON!"
    fi
    if [[ $DOPFELE == $NTEST3 ]]; then
	echo "PF electrons are OFF!"
    fi
    if [[ $DOPFELE != $NTEST3 ]]; then
	echo "PF electrons are ON!"
    fi
    if [[ -n $DATA ]]; then
	echo "Running on data!"
    fi
    if [[ -n $MC ]]; then
    	echo "Running on MC"
    	echo "MC configuration:"
    	if [[ $PYTHIA8 == $NTEST3 ]]; then
    	    echo "Special PYTHIA8 setting is OFF!"
    	fi
    	if [[ $PYTHIA8 != $NTEST3 ]]; then
    	    echo "Special PYTHIA8 setting is ON!"
    	fi
    	if [[ $ISPYTHIASHOWERED == $NTEST3 ]]; then
    	    echo "Pythia Kinematics Filter is OFF!"
    	fi
    	if [[ $ISPYTHIASHOWERED != $NTEST3 ]]; then
    	    echo "Pythia Kinematics Filter is ON!"
    	fi
    	if [[ $SHERPA == $NTEST3 ]]; then
    	    echo "Special SHERPA setting is OFF!"
    	fi
    	if [[ $SHERPA != $NTEST3 ]]; then
    	    echo "Special SHERPA setting is ON!"
    	fi
    	if [[ $SUSYPAR == $NTEST3 ]]; then
    	    echo "Special Susy Parameter setting is OFF!"
    	fi
    	if [[ $SUSYPAR != $NTEST3 ]]; then
    	    echo "Special Susy Parameter setting is ON!"
    	fi
    	if [[ $MATCHALL == $NTEST3 ]]; then
    	    echo "WARNING! MC Truth Matching only for e, mu and their neutrinos!"
    	fi


    fi
    echo "Using $GLOBALTAG as Global Tag"


    if [[ -n $JSON ]] ; then
	check_file $JSON "ERROR: JSON file $JSON does not exist"
    fi

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
    echo "Creating (recursively) output directory $STORAGE_DIR"
    srm_mkdir $OUTPUT_DIR

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
