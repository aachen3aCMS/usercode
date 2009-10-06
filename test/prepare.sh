#!/bin/sh
#
# Prepare crab submission for a given dataset
#
# Usage: ./prepare.sh <datasetpath> <version> <tag> <fast|full>
#
#   Carsten Magass, January 2009, April 2009, October 2009
#

echo ""  
echo "  -- Preparing crab submission --"
echo

if [ $# -ne 4 ]
then
  echo " ERROR  "
  echo " Wrong number of arguments !"
  echo " Usage: ./prepare.sh <datasetpath> <version> <tag> <reco|pat>"
  echo
  exit
fi


if [ $4 == "reco" ]
then
  cfg=ACSkim_cfg_reco.py
elif [ $4 == "pat" ]
then
  cfg=ACSkim_cfg_pat.py
else
  echo " ERROR : Specify 'reco' or 'pat' "
  echo
  echo " Usage: ./prepare.sh <datasetpath> <version> <tag> <reco|pat>"
  echo
  exit
fi

./mysrmls.sh output/$2/$3 >& ttt
str=`cat ttt | grep "path does not exist"`
rm -f ttt
if [ "$str" ]
then
  echo " ERROR "
  echo " Output directory (/pnfs/physik.rwth-aachen.de/cms/store/user/magass/output/$2/$3) does not exist !"
  echo ""
  echo " You might want to create it via"
  echo "   srmmkdir srm://grid-srm.physik.rwth-aachen.de:8443//pnfs/physik.rwth-aachen.de/cms/store/user/magass/output/$2[/$3]"
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

geo=`grep globaltag $cfg | cut -f1 -d":" | cut -b43-`
echo " You IMPLICITLY specified the following option "
echo "  + process.GlobalTag.globaltag : " $geo
echo ""

#if echo "$1" | grep $geo >/dev/null
#then
#  echo ""
#else
#  echo
#  echo " ERROR : Discrepancy between '$1' and '$geo' in '"$cfg"' !"
#  echo
#  exit
#fi

echo " You specified the following options "
echo "  + datasetpath     : $1 "
echo "  + user_remote_dir : output/$2/$3 "

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
echo "pset = " $cfg >> $CRABFILE
echo "total_number_of_events = -1" >> $CRABFILE
echo "events_per_job = 20000" >> $CRABFILE
echo "output_file = out.root" >> $CRABFILE
echo "" >> $CRABFILE
echo "[USER]" >> $CRABFILE
echo "return_data = 0" >> $CRABFILE
echo "email=magass@cern.ch" >> $CRABFILE
echo "copy_data = 1" >> $CRABFILE
echo "storage_element = T2_DE_RWTH" >> $CRABFILE
echo "# storage_path = /pnfs/physik.rwth-aachen.de" >> $CRABFILE
echo "user_remote_dir = output/$2/$3/" >> $CRABFILE
echo "" >> $CRABFILE
echo "[GRID]" >> $CRABFILE
echo "ce_black_list = T2_ES_IFCA,T2_TW_Taiwan,T2_US_Nebraska" >> $CRABFILE
echo "ce_white_list = T2_DE_RWTH,T2_DE_DESY,T2_US_UCSD,T2_US_Wisconsin,T2_US_MIT,T2_US_Purdue" >> $CRABFILE
echo "" >> $CRABFILE


cp $cfg $DIR/.
mv $CRABFILE $DIR/crab.cfg
echo
echo "         DONE "
echo
echo " You may now proceed with "
echo "   voms-proxy-init -voms cms:/cms/dcms -valid 164:00"
echo "   cd "$DIR
echo "   crab -create"
echo "   crab -submit"
echo "   crab -status"
echo

  





