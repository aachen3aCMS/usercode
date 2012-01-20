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
                [ -j json-file -n njobs ]

OPTIONS:
  -h               This help message
  -m 		   Running on MC
  -d               Running on data
  -s dataset       Dataset path, e.g. /DoubleMu/Run2011A-May10ReReco-v1/AOD
  -t tag           Any tag that allows you to identify this sample, e.g. "data"
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
    sed -e s/@ISDATA@/$ISDATA/g -e s/@GLOBALTAG@/$GLOBALTAG/g -e s/@QSCALE_LOW@/$QSCALELOW/g -e s/@QSCALE_HIGH@/$QSCALEHIGH/g -e s/@ELEPT@/$ELEPT/g -e s/@ELEETA@/$ELEETA/g -e s/@NELE@/$NELE/g -e s/@NMUO@/$NMUO/g -e s/@MUOPTFIRST@/$MUOPTFIRST/g -e s/@MUOPTOTHER@/$MUOPTOTHER/g -e s/@MUOETA@/$MUOETA/g -e s/@PFELEPT@/$PFELEPT/g -e s/@PFELEETA@/$PFELEETA/g -e s/@NPFELE@/$NPFELE/g -e s/@PHOPT@/$PHOPT/g -e s/@PHOETA@/$PHOETA/g -e s/@NPHO@/$NPHO/g -e s/@TAUPT@/$TAUPT/g -e s/@TAUETA@/$TAUETA/g -e s/@NTAU@/$NTAU/g -e s/@PFJETPT@/$PFJETPT/g -e s/@PFJETETA@/$PFJETETA/g -e s/@NPFJET@/$NPFJET/g -e s/@ELEPT@/$ELEPT/g -e s/@ELEETA@/$ELEETA/g -e s/@NELE@/$NELE/g -e s/@CALOJETPT@/$CALOJETPT/g -e s/@CALOJETETA@/$CALOJETETA/g -e s/@NCALOJET@/$NCALOJET/g -e s/@METCALO@/$METCALO/g -e s/@METTC@/$METTC/g -e s/@METPF@/$METPF/g -e s/@PFHTC@/$PFHTC/g -e s/@HTC@/$HTC/g -e s/@TRIGGERCONTAINS@/$TRIGGERCONTAINS/g -e s/@DOTAU@/$DOTAU/g -e s/@MUOMINV@/$MUOMINV/g -e s/@MUODMINV@/$MUODMINV/g  -e s/@DOCALOJETS@/$DOCALOJETS/g -e s/@PYTHIA8@/$PYTHIA8/g -e s/@SHERPA@/$SHERPA/g -e s/@MATCHALL@/$MATCHALL/g -e s/@SUSYPAR@/$SUSYPAR/g -e s/@FASTJET@/$FASTJET/g $CMSSWCFG >> $PYFILE
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
    MUOPTFIRST=
    MUOPTOTHER=
    NMUO=
    MUOETA=
    ELEPT=
    ELEETA=
    NELE=
    PHOPT=
    PHOETA=
    NPHO=
    PFELEPT=
    PFELEETA=
    NPFELE=
    CALOJETPT=
    CALOJETETA=
    NCALOJET=
    PFJETPT=
    PFJETETA=
    NPFJET=
    TAUPT=
    TAUETA=
    NTAU=
    METCALO=
    METPF=
    METTC=
    HTC=
    PFHTC=
    TRIGGERCONTAINS=
    MUOMINV=
    MUODMINV=
    DOTAU=
    DOCALOJETS=
    FASTJET=
    PYTHIA8=
    MATCHALL=
    SHERPA=
    SUSYPAR=
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
    SEG=`grep "muoptfirst" $CUTCFG`
    SEG2=( $SEG )
    MUOPTFIRST=${SEG2[1]}

    SEG=`grep "muoptfirst" $CUTCFG`
    SEG2=( $SEG )
    MUOPTFIRST=${SEG2[1]}

    SEG=`grep "muoptother" $CUTCFG`
    SEG2=( $SEG )
    MUOPTOTHER=${SEG2[1]}

    SEG=`grep "nmuo" $CUTCFG`
    SEG2=( $SEG )
    NMUO=${SEG2[1]}

    SEG=`grep "muoeta" $CUTCFG`
    SEG2=( $SEG )
    MUOETA=${SEG2[1]}

    SEG=`grep "elept" $CUTCFG`
    SEG2=( $SEG )
    ELEPT=${SEG2[1]}

    SEG=`grep "nele" $CUTCFG`
    SEG2=( $SEG )
    NELE=${SEG2[1]}

    SEG=`grep "eleeta" $CUTCFG`
    SEG2=( $SEG )
    ELEETA=${SEG2[1]}


    SEG=`grep "pfelept" $CUTCFG`
    SEG2=( $SEG )
    PFELEPT=${SEG2[1]}

    SEG=`grep "npfele" $CUTCFG`
    SEG2=( $SEG )
    NPFELE=${SEG2[1]}

    SEG=`grep "pfeleeta" $CUTCFG`
    SEG2=( $SEG )
    PFELEETA=${SEG2[1]}


    SEG=`grep "phopt" $CUTCFG`
    SEG2=( $SEG )
    PHOPT=${SEG2[1]}

    SEG=`grep "npho" $CUTCFG`
    SEG2=( $SEG )
    NPHO=${SEG2[1]}

    SEG=`grep "phoeta" $CUTCFG`
    SEG2=( $SEG )
    PHOETA=${SEG2[1]}


    SEG=`grep "pfelept" $CUTCFG`
    SEG2=( $SEG )
    PFELEPT=${SEG2[1]}

    SEG=`grep "npfele" $CUTCFG`
    SEG2=( $SEG )
    NPFELE=${SEG2[1]}

    SEG=`grep "pfeleeta" $CUTCFG`
    SEG2=( $SEG )
    PFELEETA=${SEG2[1]}


    SEG=`grep "calojetpt" $CUTCFG`
    SEG2=( $SEG )
    CALOJETPT=${SEG2[1]}

    SEG=`grep "ncalojet" $CUTCFG`
    SEG2=( $SEG )
    NCALOJET=${SEG2[1]}

    SEG=`grep "calojeteta" $CUTCFG`
    SEG2=( $SEG )
    CALOJETETA=${SEG2[1]}


    SEG=`grep "pfjetpt" $CUTCFG`
    SEG2=( $SEG )
    PFJETPT=${SEG2[1]}

    SEG=`grep "npfjet" $CUTCFG`
    SEG2=( $SEG )
    NPFJET=${SEG2[1]}

    SEG=`grep "pfjeteta" $CUTCFG`
    SEG2=( $SEG )
    PFJETETA=${SEG2[1]}


    SEG=`grep "taupt" $CUTCFG`
    SEG2=( $SEG )
    TAUPT=${SEG2[1]}

    SEG=`grep "ntau" $CUTCFG`
    SEG2=( $SEG )
    NTAU=${SEG2[1]}

    SEG=`grep "taueta" $CUTCFG`
    SEG2=( $SEG )
    TAUETA=${SEG2[1]}


    SEG=`grep "metcalo" $CUTCFG`
    SEG2=( $SEG )
    METCALO=${SEG2[1]}

    SEG=`grep "metpf" $CUTCFG`
    SEG2=( $SEG )
    METPF=${SEG2[1]}

    SEG=`grep "metc" $CUTCFG`
    SEG2=( $SEG )
    METTC=${SEG2[1]}


    SEG=`grep "Htc" $CUTCFG`
    SEG2=( $SEG )
    HTC=${SEG2[1]}

    SEG=`grep "PFhtc" $CUTCFG`
    SEG2=( $SEG )
    PFHTC=${SEG2[1]}

    SEG=`grep "triggerContains" $CUTCFG`
    SEG2=( $SEG )
    TRIGGERCONTAINS=${SEG2[1]}


    SEG=`grep "muoMinv" $CUTCFG`
    SEG2=( $SEG )
    MUODMINV=${SEG2[1]}

    SEG=`grep "muoDMinv" $CUTCFG`
    SEG2=( $SEG )
    MUOMINV=${SEG2[1]}

    SEG=`grep "DoTau" $CUTCFG`
    SEG2=( $SEG )
    DOTAU=${SEG2[1]}

    if [[ $DATA == $NTEST4 ]]; then
    	SEG=`grep "GlobalTagData" $CUTCFG`
    	SEG2=( $SEG )
    	GLOBALTAG=${SEG2[1]}
    fi

    if [[ $MC == $NTEST4 ]]; then
    	SEG=`grep "GlobalTagMC" $CUTCFG`
    	SEG2=( $SEG )
    	GLOBALTAG=${SEG2[1]}
    fi

    SEG=`grep "DoCaloJets" $CUTCFG`
    SEG2=( $SEG )
    DOCALOJETS=${SEG2[1]}

    SEG=`grep "SusyPar" $CUTCFG`
    SEG2=( $SEG )
    SUSYPAR=${SEG2[1]}

    SEG=`grep "FastJets" $CUTCFG`
    SEG2=( $SEG )
    FASTJET=${SEG2[1]}

    SEG=`grep "MatchAll" $CUTCFG`
    SEG2=( $SEG )
    MATCHALL=${SEG2[1]}

    SEG=`grep "Pythia8" $CUTCFG`
    SEG2=( $SEG )
    PYTHIA8=${SEG2[1]}

    SEG=`grep "Sherpa" $CUTCFG`
    SEG2=( $SEG )
    SHERPA=${SEG2[1]}

    SEG=`grep "user" $CUTCFG`
    SEG2=( $SEG )
    USER=${SEG2[1]}

    SEG=`grep "version" $CUTCFG`
    SEG2=( $SEG )
    VERSION=${SEG2[1]}

    echo "Cuts read from file"

    # sanity checks on skimming cuts. Gives warnings, does not prevent you from doing something you really want to.

    # sanity checks on skimming cuts
    if  [ $NELE = $NTEST2 ] && [ $NMUO = $NTEST2 ] && [ $NPFELE = $NTEST2 ] && [ $NPHO = $NTEST2 ] && [ $NTAU = $NTEST2 ] && [ $NPFJET = $NTEST2 ] && [ $NCALOJET = $NTEST2 ] && [ $METCALO = $NTEST ] && [ $METPF = $NTEST ] && [ $METTC = $NTEST ] ; then
    #~ if [ $NELE = $NTEST2 ] ; then
	echo " WARNING! You have spedicifed no object to cut on!"
    fi

    if [ $NELE = $NTEST2 ] && [ $ELEPT != $NTEST ]; then
	echo "WARNING! elept > 0 but nele = 0. Cut has no effect!"
    fi
    if [ $NMUO = $NTEST2 ] && [ $MUOPT != $NTEST ]; then
	echo "WARNING! muopt > 0 but nmuo = 0. Cut has no effect!"
    fi
    if [ $NPFELE = $NTEST2 ] && [ $PFELEPT != $NTEST ]; then
	echo "WARNING! pfelept > 0 but npfele = 0. Cut has no effect!"
    fi
    if [ $NTAU = $NTEST2 ] && [ $TAUPT != $NTEST ]; then
	echo "WARNING! taupt > 0 but ntau = 0. Cut has no effect!"
    fi
    if [ $NPHO = $NTEST2 ] && [ $PHOPT != $NTEST ]; then
	echo "WARNING! phopt > 0 but npho = 0. Cut has no effect!"
    fi
    if [ $NPFJET = $NTEST2 ] && [ $PFJETPT != $NTEST ]; then
	echo "WARNING! pfjetpt > 0 but npfjet = 0. Cut has no effect!"
    fi
    if [ $NCALOJET = $NTEST2 ] && [ $CALOJETPT != $NTEST ]; then
	echo "WARNING! calojetpt > 0 but ncalojet = 0. Cut has no effect!"
    fi

    if [ $NELE != $NTEST2 ]; then
	echo "Preparing Skim requiring $NELE electron(s) with pt > $ELEPT GeV and eta < $ELEETA"
    fi
    if [ $NMUO != $NTEST2 ]; then
	echo "Preparing Skim requiring $NMUO muon(s) with pt > $MUOPTFIRST GeV for the first and $MUOPTOTHER GeV for the others and eta < $MUOETA"
    fi
    if [ $NPFELE != $NTEST2 ]; then
	echo "Preparing Skim requiring $NPFELE PFelectron(s) with pt > $PFELEPT GeV and eta < $PFELEETA"
    fi
    if [ $NPHO != $NTEST2 ]; then
	echo "Preparing Skim requiring $NPHO photon(s) with pt > $PHOPT GeV and eta < $PHOETA"
    fi
    if [ $NTAU != $NTEST2 ]; then
	echo "Preparing Skim requiring $NTAU tau(s) with pt > $TAUPT GeV and eta < $TAUETA"
    fi
    if [ $NCALOJET != $NTEST2 ]; then
	echo "Preparing Skim requiring $NCALOJET calojet(s) with pt > $CALOJETPT GeV and eta < $CALOJETETA"
    fi
    if [ $NPFJET != $NTEST2 ]; then
	echo "Preparing Skim requiring $NPFJET pfjet(s) with pt > $PFJETPT GeV and eta < $PFJETETA"
    fi
    if [ $METPF != $NTEST ]; then
	echo "Preparing Skim requiring PFMET > $METPF GeV"
    fi
    if [ $METCALO != $NTEST ]; then
	echo "Preparing Skim requiring CALOMET > $METCALO GeV"
    fi
    if [ $METTC != $NTEST ]; then
	echo "Preparing Skim requiring TCMET > $METTC GeV"
    fi

    echo "Preparing crab task with following options:"
    NTEST3="False"

    if [ $DOTAU = $NTEST3 ]; then
	echo "Taus are OFF!"
    fi
    if [ $DOTAU != $NTEST3 ]; then
	echo "Taus are ON!"
    fi
    if [ $DOCALOJETS = $NTEST3 ]; then
	echo "Calojets are OFF!"
    fi
    if [ $DOCALOJETS != $NTEST3 ]; then
	echo "Calojets are ON!"
    fi
    if [[ -n $DATA ]]; then
	echo "Running on data!"
    fi
    if [[ -n $MC ]]; then
    	echo "Running on MC"
    	echo "MC configuration:"
    	if [ $PYTHIA8 = $NTEST3 ]; then
    	    echo "Special PYTHIA8 setting is OFF!"
    	fi
    	if [ $PYTHIA8 != $NTEST3 ]; then
    	    echo "Special PYTHIA8 setting is ON!"
    	fi
    	if [ $SHERPA = $NTEST3 ]; then
    	    echo "Special SHERPA setting is OFF!"
    	fi
    	if [ $SHERPA != $NTEST3 ]; then
    	    echo "Special SHERPA setting is ON!"
    	fi
    	if [ $SUSYPAR = $NTEST3 ]; then
    	    echo "Special Susy Parameter setting is OFF!"
    	fi
    	if [ $SUSYPAR != $NTEST3 ]; then
    	    echo "Special Susy Parameter setting is ON!"
    	fi
    	if [ $MATCHALL = $NTEST3 ]; then
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
