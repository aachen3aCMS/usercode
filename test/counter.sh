#!/bin/sh
#
# Count jobs, number of processed & written events for a given joboutput
#
# Usage: ./counter.sh <user> <version> <tag> 
#
#   Carsten Magass, January 2009
#                   April 2010 (update)
#

if [ $# -ne 3 ]
then
  echo
  echo " ERROR "
  echo " Wrong number of arguments !"
  echo " Usage: ./counter.sh <user> <version> <tag> "
  echo
  exit
fi

echo ""  
echo "  -- Counting Events using logfiles --"
echo


DIR=`echo "CRAB-"$3"-$2"`
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
  m=0
  TOTAL_EVENTS=0
  PASS_EVENTS=0
  OUT_EVENTS=0
  for file in $LIST
  do
#    echo $file
    SEG=`grep "segmentation" $file`
    if [ "$SEG" ]
    then
      echo "  -> segmentation violation in " $file
      let k=k+1
      let m=m+1
      continue
    fi
    EX=`grep "CMSException" $file`
    if [ "$EX" ]
    then
      echo "  -> CMSException in " $file
      let k=k+1
      let m=m+1
      continue
    fi
    file2=`echo $file | sed s/"err"/"out"/g`

    EX2=`grep "StageOutExitStatus = 60303" $file2`
    if  [ "$EX2" ]
    then
      echo "  -> File already exists in SE in " $file2
#      let k=k+1
       let m=m+1
#      continue
    fi
    SE=`grep "StageOutExitStatus = 60307" $file2`
    if [ "$SE" ]
    then
      echo "  -> Stage out problem in " $file2
      let k=k+1
      let m=m+1
      continue
    fi
    EXB=`grep "JOB_EXIT_STATUS = -1" $file2`
    if [ "$EXB" ]
    then
      echo "  -> Exit status -1 in " $file2
      let k=k+1
      let m=m+1
      continue
    fi
    EXC=`grep "StageOutExitStatus = 8001" $file2 | grep -v "record"`
    if [ "$EXC" ]
    then
      echo "  -> Exit status 8001 in " $file2
      let k=k+1
      let m=m+1
      continue
    fi
    EXA=`grep "EXECUTABLE_EXIT_STATUS = 8001" $file2`

    if [ "$EXA" ]
    then
      echo "  -> Exit status 8001 in " $file2
      let k=k+1
      let m=m+1
      continue
    fi	
    EXF=`grep "EXECUTABLE_EXIT_STATUS = 8020" $file2`

    if [ "$EXF" ]
    then
      echo "  -> Exit status 8020 in " $file2
      let k=k+1
      let m=m+1
      continue
    fi	

    EVENTS1=`grep "Events total" $file2 | awk '{print $5}'`
    EVENTS2=`grep "of events" $file2 | awk '{print $6}'`
    EVENTS3=`grep "of events" $file2 | awk '{print $10}'`
#    echo $file "  "$EVENTS1"|"$EVENTS2"|"$EVENTS3
    TOTAL_EVENTS=$(echo "$TOTAL_EVENTS+$EVENTS1" | bc)
    PASS_EVENTS=$(echo "$PASS_EVENTS+$EVENTS2" | bc)
    OUT_EVENTS=$(echo "$OUT_EVENTS+$EVENTS3" | bc)
    let k=k+1
    let j=j+1
  done
  echo
  echo '  checked '$j' (out of '$k') CMSSW_*.stderr files '
  echo '       -->  '$m' files contain error messages '
  echo '       ==>  #events read    : ' $TOTAL_EVENTS
  echo '       ==>  #events presel  : ' $PASS_EVENTS
  echo '       ==>  #events written : ' $OUT_EVENTS
  echo
else
  echo "  ERROR "
  echo "  Directory ($DIR) does not contain log files of type CMSSW_*.stderr."
  echo
fi

  LAST=`ls -td $DIR/crab_* | head -n1`
  echo '  compare with crab submit '$LAST ': '
#  JOBS=`grep "can" $1/crab_*/log/crab.log | awk '{print $1}'`
#  EV=`grep "can" $1/crab_*/log/crab.log | awk '{print $6}'`
if [ -f $LAST/log/crab.log ]
then
  JOBS=`grep "job(s)" $LAST/log/crab.log | awk '{print $4}'`
  EV=`grep "job(s)" $LAST/log/crab.log | awk '{print $9}'`
  echo '     ' $JOBS 'jobs over' $EV 'events'
else
  echo
  echo "  ERROR"
  echo "  Log file ($LAST/log/crab.log) does not exist ! " 
fi

echo ""  
echo "  -- Counting Files in dcache --"
echo
f=`echo "/pnfs/physik.rwth-aachen.de/cms/store/user/"$1"/output/"$2"/"$3`

#echo $f
#n=`ls $f/*root | wc -l`
#n=`srmls srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/magass/output/$1/$2 | grep -v WARNING | grep pnfs | grep -v SURL | grep root | wc -l`

n=`./mysrmls.sh $1 output/$2/$3 | grep root | wc -l`
echo "  found "$n" files in $f"
echo
nu=`./mysrmls.sh $1 output/$2/$3 | grep root | awk -F_ '{print $1 "_" $2 }' | sort | uniq -d | sort -u | wc -l`
if [ $nu -gt 0 ]
then
  echo "     -> the following files appear more than once:"
  list=`./mysrmls.sh $1 output/$2/$3 | grep root | awk -F_ '{print $1 "_" $2 }' | sort | uniq -d | sort -u`
  for s in $list
  do
    echo "        " $s
  done
else
  echo "     -> no duplicate files found :-)"
fi

echo



