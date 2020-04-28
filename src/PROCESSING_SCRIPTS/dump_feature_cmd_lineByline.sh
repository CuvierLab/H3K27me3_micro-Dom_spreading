#!/bin/bash

# source conf file
while getopts "c:" option
do
case $option in
    c)
        conffile=$OPTARG
        ;;
    \?)
        echo "-c emplacement of configuration file"
        exit 1
        ;;
    h)
    echo "-c emplacement of configuration file"
    ;;
esac
done


. ${conffile}

echo ${myhic}
echo ${dumpdir}
# Declare constant
declare -i cnt=0
declare -i max=4


mkdir -p ${dumpdir}
for file in $(ls $mydir)
do
  if (( "$cnt" <= "$max" ))
  then
    myfile=$mydir/$file
    echo $myfile
    ./dump_features_lineByline.sh -f $myfile -h ${myhic} -o $dumpdir/ -r ${res} &
    ((cnt++))
  else
    wait
    myfile=$mydir/$file
    ./dump_features_lineByline.sh -f $myfile -h ${myhic} -o $dumpdir/ -r ${res} &
    declare -i cnt=0
  fi
done

echo "end"
