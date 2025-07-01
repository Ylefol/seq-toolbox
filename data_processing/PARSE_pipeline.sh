#!/bin/usr/env bash

#Define Variables
################################################################################
experiment=$1
fastqDir=$2
nthreads=$3
genome_path=$4

#May have to be automated
script_dir="/home/yohanl/A_Projects/seq-toolbox/data_processing/"

#Set timer
start=`date +%s`

#Set output directories
################################################################################
datadir=./results && mkdir -p $datadir/$experiment/



#Find R1 and R2
################################################################################
find $fastqDir -maxdepth 1 -type f -iname "${experiment}*.fastq.gz" | sort > filelist_"$experiment"

numb_of_files=`cat filelist_"$experiment" | wc -l`
r1=`cat filelist_"$experiment" | head -n1`
r2="" #Create and empty R2, will be overwritten if a real one exists
echo "R1 is " $r1
if [ $numb_of_files == 2 ]; then
    r2=`cat filelist_"$experiment" | tail -n1`
    echo "R2 is " $r2
fi
if [ $numb_of_files == 0 ]; then
    echo "R1 does not exist "
    exit 1
fi
rm filelist_"$experiment"

# Run split-pipe
################################################################################
PARSE_output=$datadir/$experiment/
split-pipe --mode all --chemistry v3 --kit WT --nthreads $nthreads --fq1 $r1 --fq2 $r2 --output_dir $PARSE_output --genome_dir $genome_path --sample all-well A1-D12


end=`date +%s`
runtime=$((end-start))
echo "Pipeline completed for $experiment in seconds for combination"
