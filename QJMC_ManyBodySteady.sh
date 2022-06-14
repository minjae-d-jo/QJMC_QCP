#!/bin/bash

number_of_time_step=$1
N=$2
omega=$3
nRepeat=$4
uid=$5

[[ $# -ne 5 ]] && echo "Error: You need to give 5 arguments." && exit 1


params=${number_of_time_step}_${N}_${omega}_${nRepeat}_${uid}
logFile=Log/QJMC_$params.table
outputFile=Result/QJMC_$params.table

exec 2>> $logFile

echo "$(date +%Y-%m-%d_%H:%M:%S): start $$ $(hostname)" >> $logFile

Bin/QJMC_ManyBodySteady_${N} $number_of_time_step $N $omega $nRepeat >> $outputFile #2>> $outputFile2
ret=$?;
echo -n "$(date +%Y-%m-%d_%H:%M:%S): " >> $logFile
if [ $ret -eq 0 ]; then
	echo "normal exit"
else
	echo "abnormal exit: $ret"
fi >> $logFile
