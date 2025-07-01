#!/usr/bin/env bash

r1=$1
r2=$2
experiment=$3
out=$4
refgen=$5
numcores=$6


echo "################################################################################"
echo 'Aligning reads to the reference genome ...'
echo "################################################################################"
if [ "$r2" = "NONE" ]; then #SE
	echo "Mapping SE with bowtie2"
	bowtie2 -p $numcores -x $refgen -U $r1 -S "$out/$experiment.sam"
else
	echo "Mapping PE with bowtie2"
	bowtie2 -p $numcores -x $refgen -1 $r1 -2 $r2 -S "$out/$experiment.sam"
fi


echo "################################################################################"
echo "Sorting"
echo "################################################################################"
samtools view -bhS $out/$experiment".sam" | samtools sort -O bam -@ $numbproc - -o $out/$experiment"_sorted.bam" 

