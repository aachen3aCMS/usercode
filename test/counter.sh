#!/bin/sh
#
# Count jobs, number of processed & written events for a given joboutput
#
# Usage: ./counter.sh <version> <tag> 
#
#   Carsten Magass, January 2009
#

if [ $# -ne 2 ]
then
  echo
  echo " ERROR "
  echo " Wrong number of arguments !"
  echo " Usage: ./counter.sh <version> <tag> "
  echo
  exit
fi

echo ""  
echo "  -- Counting Files in dcache --"
echo
f=`echo "/pnfs/physik.rwth-aachen.de/cms/store/user/magass/output/"$1"/"$2`
#n=`ls $f/*root | wc -l`
#n=`srmls srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/magass/output/$1/$2 | grep -v WARNING | grep pnfs | grep -v SURL | grep root | wc -l`
n=`mysrmls.sh output/$1/$2 | grep root | wc -l`
echo "  found "$n" files in $f"
echo

echo ""  
echo "  -- Counting Events using logfiles --"
echo


DIR=`echo "CRAB-"$2"-$1"`
if [ ! -d $DIR ]
then
  echo " ERROR "
  echo " Directory with logfiles ($DIR) does not exist !"
  echo ""
  exit
fi

DIR3=`ls -td $DIR/crab_*/res | head -n1`
if [ ! -d $DIR3 ]
then
  echo " ERROR "
  echo " Directory with logfiles ($DIR3) does not exist !"
  echo ""
  exit
fi


LIST=`find $DIR/crab_*/res -type f -print0 | xargs -0 ls | grep stderr`

if [ "$LIST" ]
then
  j=0
  k=0
  IN_EVENTS=0
  OUT_EVENTS=0
  for file in $LIST
  do
#    echo $file
    SEG=`grep "segmentation" $file`
    if [ "$SEG" ]
    then
      echo "  -> segmentation violation in " $file
      let k=k+1
      continue
    fi
    EX=`grep "CMSException" $file`
    if [ "$EX" ]
    then
      echo "  -> CMSException in " $file
      let k=k+1
      continue
    fi
    file2=`echo $file | sed s/"err"/"out"/g`
    EX2=`grep "60303" $file2`
    if  [ "$EX2" ]
    then
      echo "  -> File already exists in SE in " $file2
      let k=k+1
#      continue
    fi
    SE=`grep "60307" $file2`
    if [ "$SE" ]
    then
      echo "  -> Stage out problem in " $file2
      let k=k+1
      continue
    fi
    EXB=`grep "JOB_EXIT_STATUS = -1" $file2`
    if [ "$EXB" ]
    then
      echo "  -> Exit status -1 in " $file2
      let k=k+1
      continue
    fi
    EXC=`grep "8001" $file2`
    if [ "$EXC" ]
    then
      echo "  -> Exit status -1 in " $file2
      let k=k+1
      continue
    fi

    EVENTS1=`grep "events" $file | awk '{print $6}'`
    EVENTS2=`grep "events" $file | awk '{print $10}'`
#    echo $file "  " $EVENTS1 "  " $EVENTS2
    IN_EVENTS=$(echo "$IN_EVENTS+$EVENTS1" | bc)
    OUT_EVENTS=$(echo "$OUT_EVENTS+$EVENTS2" | bc)
    let k=k+1
    let j=j+1
  done
  echo
  echo '  checked '$j' (out of '$k') CMSSW_*.stderr files '
  echo '       ==>  #events read    : ' $IN_EVENTS
  echo '       ==>  #events written : ' $OUT_EVENTS
  echo
  LAST=`ls -td $DIR/crab_* | head -n1`
  echo '  compare with crab submit '$LAST ': '
#  JOBS=`grep "can" $1/crab_*/log/crab.log | awk '{print $1}'`
#  EV=`grep "can" $1/crab_*/log/crab.log | awk '{print $6}'`
  JOBS=`grep "job(s)" $LAST/log/crab.log | awk '{print $4}'`
  EV=`grep "job(s)" $LAST/log/crab.log | awk '{print $9}'`
  echo '     ' $JOBS 'jobs over' $EV 'events'
echo
else
  echo " Directory ($DIR) does not contain log files of type CMSSW_*.stderr."
  echo
fi
