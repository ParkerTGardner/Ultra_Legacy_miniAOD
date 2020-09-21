#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters"
    echo "<Inputs File> <nJobs>"
    exit
fi

inDir=$1
nJobs=$2

for i in $(seq 0 $(($nJobs-1)) );
do
    ./bin/MyClass.exe $inDir $i $nJobs > 'unmergedOutputs/job_MyClass_'$i'.log' 2>&1 & 
done
