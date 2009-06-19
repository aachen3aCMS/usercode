if [ $# -eq 1 ]
then
  offset=0
elif [ $# -eq 2 ]
then
  offset=$2
else
  exit
fi

s=`echo $1`
    #  echo $s
if [ $s == '.' ]
then
#  echo "ls /pnfs/physik.rwth-aachen.de/cms/store/user/magass/$1"
    srmls -offset=$offset srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/magass/$1 | grep -v WARNING | grep pnfs | awk '{print $2}' | sed s/"\.\/"/""/g | cut -c51- | grep "/" 
  else
    srmls -offset=$offset srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/magass/$1 | grep -v WARNING | grep pnfs | awk '{print $2}' | sed s/"\.\/"/""/g | cut -c51- | grep "/" | sed s:$s:: | cut -c2-
  fi
  count=`srmls -offset=$offset srm://grid-srm.physik.rwth-aachen.de:8443/pnfs/physik.rwth-aachen.de/cms/store/user/magass/$1 | grep -v WARNING | grep pnfs | grep -v SURL | sed s:/pnfs/physik.rwth-aachen.de/cms/store/user/magass/::g  | grep root | wc -l`
#  echo "  #Files : " $count
  if [ $count -eq 999 ]
  then
    let offset=offset+999
    mysrmls.sh $1 $offset
fi 

