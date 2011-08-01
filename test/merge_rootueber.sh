#!/bin/sh
#
# Merge multiple root files located in dcache
#
# Usage: ./merge_root.sh <username> <version> <sample>
#
#   Carsten Magass, April 2009
#                   April 2010 (update)
#   Klaas Padeken Aug 2011 (upadte)
#

if [ $# -lt 3 ] || [ $# -gt 4 ] 
then
  echo "Usage: ./merge_root.sh <username> <version> <sample> [number of files merged]"
  exit 
fi

COUNTER=0
id=1
final=0
num=10

if [ $# -eq 4 ]
then 
  num=$4
fi

if [ -f templist ]
then
  rm -f templist
fi

PP=`pwd`
tmp=`./mysrmls.sh $1 output/$2/$3 | grep root > templist`
nfiles=`wc -l templist | awk '{print $1}'`

echo ""
echo "  COPY & MERGE root files from dcache"
echo "  -----------------------------------"
echo ""
echo "    Working directory   : $PP "
echo "    Temporary file list : templist "
echo "    Number of files     : $nfiles "
echo "    Number of subsets   : $[$[nfiles-1]/$num+1] " 
echo ""

for file in $( cat templist )
do
  final=1
  s=`echo "$s dcap://grid-se110.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$file "`
  s1=`echo $s1 $file`
  
  #echo "srmcp srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$file  file:///$PP/$file"
  #srmcp srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/$file  file:///$PP/$file > /dev/null
  echo "uberftp grid-ftp.physik.rwth-aachen.de active; cd /pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/; get $file"
  uberftp grid-ftp.physik.rwth-aachen.de  "active; cd /pnfs/physik.rwth-aachen.de/cms/store/user/$1/output/$2/$3/; get $file"> /dev/null
  let COUNTER=COUNTER+1
  rest=$[$COUNTER%$num]
  if [ $rest -eq 0 ]
  then
    out=`echo $3'_'$id'.root'`
    echo 
    echo "hadd -f7 "$out 
    echo $s1
    hadd -f7 $out $s1 
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
  echo "hadd -f7 "$out 
  echo $s1   
  echo
  hadd -f7 $out $s1
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
rm -f templist

#rm -f $3_*



