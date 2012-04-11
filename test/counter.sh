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


LIST=`find $DIR/crab_*/res -type f -print0 | xargs -0 ls | grep stdout`

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
#    echo "Checking $file..."
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
    EEX=`grep "StageOutExitStatus = 60317" $file`
    if [ "$EEX" ]
    then
      echo "  -> 60317 - Forced timeout for stuck stage out for file " $file
      let k=k+1
      let m=m+1
      continue
    fi

    EX2=`grep "StageOutExitStatus = 60303" $file`
    if  [ "$EX2" ]
    then
      echo "  -> 60303 - File already exists in SE in " $file
#      let k=k+1
       let m=m+1
#      continue
    fi
    SE=`grep "StageOutExitStatus = 60307" $file`
    if [ "$SE" ]
    then
      echo "  -> 60307 - Stage out problem in " $file
      let k=k+1
      let m=m+1
      continue
    fi
    EXB=`grep "JOB_EXIT_STATUS = -1" $file`
    if [ "$EXB" ]
    then
      echo "  -> Exit status -1 in " $file
      let k=k+1
      let m=m+1
      continue
    fi
    EXC=`grep "StageOutExitStatus = 8001" $file | grep -v "record"`
    if [ "$EXC" ]
    then
      echo "  -> Exit status 8001 in " $file
      let k=k+1
      let m=m+1
      continue
    fi
    EXA=`grep "EXECUTABLE_EXIT_STATUS = 8001" $file`

    if [ "$EXA" ]
    then
      echo "  -> Exit status 8001 in " $file
      let k=k+1
      let m=m+1
      continue
    fi	
    EXF=`grep "EXECUTABLE_EXIT_STATUS = 8020" $file`

    if [ "$EXF" ]
    then
      echo "  -> Exit status 8020 in " $file
      let k=k+1
      let m=m+1
      continue
    fi
    # Check requirement: SUSYana must have finished
    FIN=`grep "SusyACSkimAnalysis:ACSkimAnalysis@endJob" $file`
    if [ "$FIN" == "" ] ; then
	echo "Endjob not found in file $file"
      let k=k+1
      let m=m+1
      continue
    fi

    EVENTS1=`grep "Events total" $file | awk '{print $5}'`
    EVENTS2=`grep "of events" $file | awk '{print $6}'`
    EVENTS3=`grep "of events" $file | awk '{print $10}'`
#    echo $file "  "$EVENTS1"|"$EVENTS2"|"$EVENTS3
    let TOTAL_EVENTS=$TOTAL_EVENTS+$EVENTS1
    let PASS_EVENTS=$PASS_EVENTS+$EVENTS2
    let OUT_EVENTS=$OUT_EVENTS+$EVENTS3
    let k=k+1
    let j=j+1
  done
  echo
  echo '  checked '$j' (out of '$k') CMSSW_*.stdout files '
  echo '       -->  '$m' files contain error messages '
  echo '       ==>  #events read    : ' $TOTAL_EVENTS
  echo '       ==>  #events presel  : ' $PASS_EVENTS
  echo '       ==>  #events written : ' $OUT_EVENTS
  echo
else
  echo "  ERROR "
  echo "  Directory ($DIR) does not contain log files of type CMSSW_*.stdout."
  echo
fi

  LAST=`ls -td $DIR/crab_* | head -n1`
  echo '  compare with crab submit '$LAST ': '
#  JOBS=`grep "can" $1/crab_*/log/crab.log | awk '{print $1}'`
#  EV=`grep "can" $1/crab_*/log/crab.log | awk '{print $6}'`
if [ -f $LAST/log/crab.log ]
then
  JOBS=`grep "Total of" $LAST/log/crab.log | grep submitted | awk '{print $6}' | head -n1`
  EV=`grep "(s)" $LAST/log/crab.log | awk '{print $9}'`
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
  echo "     -> the following files appear more than once and have to be deleted using srmrm:"
  echo
  list=`./mysrmls.sh $1 output/$2/$3 | grep root | awk -F_ '{print $1 "_" $2 }' | sort | uniq -d | sort -u`
  for s in $list
  do
    file=`./mysrmls.sh $1 output/$2/$3 | grep $s'_'`
    echo "        " $s : $file
    good=`grep $s'_' $DIR/crab_*/res/*out | grep newL  | awk '{print $4}' | awk -F'out_' '{print "out_"$2}'`
    for ff in $file
    do
      if [ "$ff" != "$good" ]
      then
        echo "            srmrm srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$ff "
        echo
      fi
    done

#    file=`./mysrmls.sh $1 output/$2/$3 | grep $s'_'`
#    ffile=`grep $s CRAB-zz-v62/crab_*/res/*out | grep newL  | awk '{print $4}' | awk -F'out_' '{print "out_"$2}'`
#    echo $file " => " $ffile
  done
else
  echo "     -> no duplicate files found :-)"
fi

echo



