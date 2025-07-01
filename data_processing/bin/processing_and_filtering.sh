#!/usr/bin/env bash

out=$1
aux=$2
experiment=$3
numcores=$4

echo "################################################################################"
echo "indexing"
echo "################################################################################"
samtools index -@ $numcores $out/$experiment"_sorted.bam" $out/$experiment"_sorted.bam.bai" 

echo "################################################################################"
echo "Running collate"
echo "################################################################################"
#Shuffles and groups reads together by their names
samtools collate -@ $numcores -o $aux/$experiment"_collate.bam" $out/$experiment"_sorted.bam"

echo "################################################################################"
echo "Running fixmate"
echo "################################################################################"
#a tool that can fill in information (insert size, cigar, mapq) about paired end reads onto the corresponding other read
samtools fixmate -@ $numcores -m $aux/$experiment"_collate.bam" $aux/$experiment"_fixmate.bam"

echo "################################################################################"
echo "Running position Sort"
echo "################################################################################"
#Sorts based on position
samtools sort -@ $numcores -o $aux/$experiment"_positionsort.bam" $aux/$experiment"_fixmate.bam"

echo "################################################################################"
echo "Running mardup"
echo "################################################################################"
#Mark/remove duplicates - We need to run collate, fixmate, and sort to be able to run markdup correctly
#The duplicates are removed thanks to the -r parameter
samtools markdup -@ $numcores -r -s $aux/$experiment"_positionsort.bam" $out/$experiment"_dedup.bam"
samtools index $out/$experiment"_dedup.bam" $out/$experiment"_dedup.bam.bai"
samtools flagstat $out/$experiment"_sorted.bam" > $out/$experiment".flagstat"

echo "Done filtering"
