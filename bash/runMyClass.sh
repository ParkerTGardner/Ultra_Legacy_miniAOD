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
    #./bin/MyClass.exe $inDir $i $nJobs > 'src/testNew/job_MyClass_'$i'.log' 2>&1 & 
    #./bin/MyClass.exe $inDir $i $nJobs > 'src/Full2018/job_MyClass_'$i'.log' 2>&1 & 
    #./bin/MyClass.exe $inDir $i $nJobs > 'src/Full2017/job_MyClass_'$i'.log' 2>&1 & 
    ./bin/MyClass.exe $inDir $i $nJobs > 'src/Full2016/job_MyClass_'$i'.log' 2>&1 & 
done
