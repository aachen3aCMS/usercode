#!/bin/sh

if [ $# -ne 3 ]
then
  echo
  echo " Usage: ./create_root_list.sh <version> <sample> <filelist>"
  echo
  exit 
fi

list=`grep out_ CRAB-$2-$1/crab_*/res/*out | grep LFN | grep store | awk '{print$2}' `

if [ -f $3 ]
then
  echo
  echo " ERROR : File list $3 already exists !"
  echo
  exit
else
  touch $3
fi
echo
echo " ... creating list named $3"

for file in $list
do
#  echo $file	
  cind=`echo $file | awk '{print index($0,"out_")}'`	
  newf=`echo $file | cut -b $cind-`
  echo $newf >> $3
done

nfiles=`cat $3 | wc -l`
echo "   -> contains $nfiles files    "
echo
