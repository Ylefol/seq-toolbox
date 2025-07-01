#!/usr/bin/env bash

out=$1
experiment=$2
genome=$3
numbproc=$4

echo "################################################################################"
echo 'Running Feature counts'
echo "################################################################################"
gtf_files=($(find "$genome" -maxdepth 1 -type f -name "*.gtf"))
if [ $numb_of_files == 2 ]; then
	echo "Feature counts (PE)"
	featureCounts -T $numbproc -p --countReadPairs -t exon -g gene_name -a $gtf_files -o $out/counts_raw.txt "$out/STAR_res/Aligned.sortedByCoord.out.bam"
else
	echo "Feature counts (SE)"
	featureCounts -T $numbproc -t exon -g gene_name -a $gtf_files -o $out/counts_raw.txt "$out/STAR_res/Aligned.sortedByCoord.out.bam"
fi


# Format the ouput
################################################################################
echo "################################################################################"
echo 'Formatting the counts ouput'
echo "################################################################################"
cut -f1,7-8 $out/counts_raw.txt > $out/$experiment.counts
sed -i '1,2d' $out/$experiment.counts
