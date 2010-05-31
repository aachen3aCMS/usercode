#!/bin/sh
#
# Merge multiple root files located in dcache
#
# Usage: ./merge_root_list.sh <username> <version> <sample>
#
#   Carsten Magass, May 2010
#

if [ $# -ne 4 ]
then
  echo
  echo " Usage: ./merge_root_list.sh <username> <version> <sample> <filelist>"
  echo
  exit 
fi

COUNTER=0
id=1
final=0

if [ -f $4 ]
then
  f=`echo 1`
else 
  echo
  echo " ERROR :  File $4 does not exit !" 
  echo
  exit
fi

PP=`pwd`
#tmp=`./mysrmls.sh $1 output/$2/$3 | grep root > templist`
nfiles=`wc -l $4 | awk '{print $1}'`

echo ""
echo "  COPY & MERGE root files from dcache"
echo "  -----------------------------------"
echo ""
echo "    Working directory   : $PP "
echo "    File list           : $4 "
echo "    Number of files     : $nfiles "
echo "    Number of subsets   : $[$nfiles/10+1] " 
echo ""

for file in $( cat $4 )
do
  final=1
  s=`echo "$s dcap://grid-se110.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$file "`
  s1=`echo $s1 $file`
  echo "srmcp srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$file  file:///$PP/$file"
  srmcp srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$file  file:///$PP/$file > /dev/null
  let COUNTER=COUNTER+1
  rest=$[$COUNTER%10]
  if [ $rest -eq 0 ]
  then
    out=`echo $3'_'$id'.root'`
    echo 
    echo "hadd "$out
    echo $s1
    hadd $out $s1
    s1=''
    out=`echo $3'_'$id'.root'`
    let id=id+1
    final=0
    rm -f out_*.root
  fi
done

if [ $final -eq '1' ]
then
  out=`echo $3'_'$id'.root'`
  echo 
  echo "hadd "$out
  echo $s1   
  echo
  hadd $out $s1
  rm -f out_*.root
fi

if [ $id -eq '1' ]
then
  out2=`echo $3'.root'`
  echo "mv $out -> $out2"
  echo
  mv $out $out2
fi
let id=id-1
#DIR=`echo $3"_$id.root"`
#if [ -e $DIR ]
#then
#  file=`echo $3".root"`
#  str=`echo $3"_*.root"`
#  echo "hadd" $file $str
#  hadd $file $str
#fi 

#rm -f $3_*



