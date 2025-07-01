#!/usr/bin/env bash


out=$1
experiment=$2
resloc=$3
effgensize=$4
numcores=$5


#Create bigwig
################################################################################
echo "################################################################################"
echo "creating bigwig"
echo "################################################################################"

echo "bigwig for sorted"
bamCoverage -p $numcores -b $out/$experiment"_sorted.bam" -o $resloc/$experiment"_sorted_RPKM.bw" --normalizeUsing RPKM --effectiveGenomeSize $effgensize --ignoreDuplicates --extendReads --ignoreForNormalization X M
echo "################################################################################"

echo "bigwig for dedup"
bamCoverage -p $numcores -b $out/$experiment"_dedup.bam" -o $resloc/$experiment"_dedup_RPKM.bw" --normalizeUsing RPKM --effectiveGenomeSize $effgensize --ignoreDuplicates --extendReads --ignoreForNormalization X M
echo "################################################################################"

#Create bed file
################################################################################
echo "################################################################################"
echo "creating bedfile"
echo "################################################################################"
bedtools bamtobed -i $out/$experiment"_dedup.bam" > $resloc/$experiment"_dedup.bed"
bedtools bamtobed -i $out/$experiment"_sorted.bam" > $resloc/$experiment"_sorted.bed"


#Run MACS3 - peak calling
################################################################################
echo "################################################################################"
echo "Peak calling for $experiment"
echo "################################################################################"

macs3 callpeak -t $out/$experiment"_dedup.bam" -n $experiment"_dedup" --outdir $resloc
macs3 callpeak -t $out/$experiment"_sorted.bam" -n $experiment"_sorted" --outdir $resloc

