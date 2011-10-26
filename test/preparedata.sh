#!/bin/sh
#
# Prepare crab submission for a given dataset
#
# Usage: ./prepare.sh <user> <datasetpath> <version> <tag> <cfg> <globaltag> [ <pthat_low> <pthat_high> ]
#
#   Carsten Magass, January 2009, April 2009, October 2009, November 2009, December 2009, April 2010
#

echo ""  
echo "  -- Preparing crab submission --"
echo

if [ $# -le 7 ] || [ $# -gt 9 ] 
then
  echo " ERROR  "
  echo " Wrong number of arguments !"
  echo " Usage: ./prepare.sh <user> <datasetpath> <version> <tag> <cfg> <globaltag> <json> <jobs>"
  echo
  exit
fi

if [ -a $5 ]
then
  dummy=1
else
  echo " ERROR : " $5 " does not exist !" 
  echo
  exit
fi 

./mysrmls.sh $1 output/$3/$4 >& ttt
str=`cat ttt | grep "does not exist"`
rm -f ttt
if [ "$str" ]
then
  echo " ERROR "
  echo " Output directory (/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$3/$4) does not exist !"
  echo ""
  echo " You might want to create it via"
  echo "   srmmkdir srm://grid-srm.physik.rwth-aachen.de:8443//pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$3/$4"
  echo
  exit
fi

TEMPDATE=`date +%Y_%m_%e-%k_%M_%S `
DIR=`echo "CRAB-"$4"-"$3`

if [ -d $DIR ]
then
  echo " ERROR "
  echo " Working directory ($DIR) already exists !"
  echo
  exit
fi

echo " You specified the following options "
echo "  + user                        : $1 "
echo "  + configuration               : $5 "
echo "  + datasetpath                 : $2 "
echo "  + user_remote_dir             : output/$3/$4 "
echo "  + process.GlobalTag.globaltag : $6"

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
echo "use_server = 1" >> $CRABFILE
echo "#server_name = cern" >> $CRABFILE
echo "" >> $CRABFILE
echo "[CMSSW]" >> $CRABFILE
echo "# dbs_url = http://cmsdbsprod.cern.ch/cms_dbs_ph_analysis_02/servlet/DBSServlet" >> $CRABFILE
echo "datasetpath = $2" >> $CRABFILE
echo "pset = " $5 >> $CRABFILE
#echo "total_number_of_events = -1" >> $CRABFILE
#echo "events_per_job = 20000" >> $CRABFILE
echo "lumi_mask = " $7 >> $CRABFILE
echo "total_number_of_lumis = -1" >> $CRABFILE
echo "number_of_jobs = "$8 >> $CRABFILE

echo "output_file = out.root" >> $CRABFILE
echo "" >> $CRABFILE
echo "[USER]" >> $CRABFILE
echo "return_data = 0" >> $CRABFILE
echo "email=$1@cern.ch" >> $CRABFILE
echo "copy_data = 1" >> $CRABFILE
echo "storage_element = T2_DE_RWTH" >> $CRABFILE
echo "user_remote_dir = output/$3/$4/" >> $CRABFILE
echo "check_user_remote_dir = 0" >> $CRABFILE
echo "" >> $CRABFILE
echo "[GRID]" >> $CRABFILE
echo "#ce_white_list = T2_DE_RWTH,T2_DE_DESY,T2_US_UCSD,T2_US_Wisconsin,T2_US_MIT,T2_US_Florida" >> $CRABFILE
echo "#to run without server: (use_server = 0)">>$CRABFILE
echo "#group = dcms" >> $CRABFILE
echo "#max_wall_clock_time = 1400" >> $CRABFILE
echo "#max_cpu_time = 1400" >> $CRABFILE
echo "#rb = CERN" >> $CRABFILE
echo "#se_black_list = T0,T1" >> $CRABFILE
echo "" >> $CRABFILE

sed s/"globaltag = cms.string('')"/"globaltag = cms.string('$6')"/g < $5 > $DIR/temp.txt

mv $DIR/temp.txt $DIR/$5
 
mv $CRABFILE $DIR/crab.cfg
cp $7 $DIR/
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

  





