#!/bin/sh
#
# Prepare crab submission for a given dataset
#
# Usage: ./prepare.sh <datasetpath> <version> <tag>
#
#   Carsten Magass, January 2009
#

echo ""  
echo "  -- Preparing crab submission --"
echo

if [ $# -ne 3 ]
then
  echo " ERROR  "
  echo " Wrong number of arguments !"
  echo " Usage: ./prepare.sh <datasetpath> <version> <tag> "
  echo
  exit
fi

if [ ! -d /pnfs/physik.rwth-aachen.de/dcms/magass/output/$2/$3 ]
then
  echo " ERROR "
  echo " Output directory (/pnfs/physik.rwth-aachen.de/dcms/magass/output/$2/$3) does not exist !"
  echo ""
  echo " You might want to create it via"
  echo "   srmmkdir srm://grid-srm.physik.rwth-aachen.de:8443//pnfs/physik.rwth-aachen.de/dcms/magass/output/$2[/$3]"
  echo
  exit
fi

TEMPDATE=`date +%Y_%m_%e-%k_%M_%S `
DIR=`echo "CRAB-"$3"-"$2`

if [ -d $DIR ]
then
  echo " ERROR "
  echo " Working directory ($DIR) already exists !"
  echo
  exit
fi

echo " You specified the following options "
echo "  + datasetpath     : $1 "
echo "  + user_remote_dir : dcms/magass/output/$2/$3 "


if [ ! -d $DIR ] 
then
  mkdir $DIR
fi

CRABFILE=tempcrab.cfg

if [ -e $CRABFILE ] 
then
  rm $CRABFILE
fi  
touch $CRABFILE

echo "[CRAB]" >> $CRABFILE
echo "jobtype = cmssw" >> $CRABFILE
echo "scheduler = glite" >> $CRABFILE
echo "server_name = cern" >> $CRABFILE
echo "" >> $CRABFILE
echo "[CMSSW]" >> $CRABFILE
echo "datasetpath = $1" >> $CRABFILE
echo "pset = ACSkim_cfg.py" >> $CRABFILE
echo "total_number_of_events = -1" >> $CRABFILE
echo "events_per_job = 10000" >> $CRABFILE
echo "output_file = out.root" >> $CRABFILE
echo "" >> $CRABFILE
echo "[USER]" >> $CRABFILE
echo "return_data = 0" >> $CRABFILE
echo "email=magass@cern.ch" >> $CRABFILE
echo "copy_data = 1" >> $CRABFILE
echo "storage_element = grid-srm.physik.rwth-aachen.de" >> $CRABFILE
echo "storage_path = /pnfs/physik.rwth-aachen.de" >> $CRABFILE
echo "user_remote_dir = dcms/magass/output/$2/$3/" >> $CRABFILE
echo "" >> $CRABFILE
echo "[EDG]" >> $CRABFILE
echo "ce_black_list = ucsd.edu" >> $CRABFILE
echo "" >> $CRABFILE


cp ACSkim_cfg.py $DIR/.
mv $CRABFILE $DIR/crab.cfg
echo
echo "         DONE "
echo
echo " You may now proceed with "
echo "   cd "$DIR
echo "   voms-proxy-init -voms cms:/cms/dcms"
echo "   crab -create"
echo "   crab -submit"
echo "   crab -status"
echo

  





