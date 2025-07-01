#!/usr/bin/env bash

r1=$1
r2=$2
out=$3
experiment=$4
genome=$5
numbproc=$6
script_dir=$7


#Create location file
STAR_loc=$script_dir/STAR_res/ && mkdir -p $STAR_loc

echo "################################################################################"
echo 'Aligning reads to the reference genome ...'
echo "################################################################################"

if [ "$r2" = "NONE" ]; then #SE

	echo "Aligning with STAR (SE) - logs will be in a seperate log file"
	STAR --runThreadN $numbproc \
	     --genomeDir "$genome" \
	     --readFilesIn "$r1" \
	     --readFilesCommand zcat \
	     --outFileNamePrefix $STAR_loc \
	     --outSAMtype BAM SortedByCoordinate \
	     --outSAMunmapped Within \
	     --outFilterType BySJout \
	     --outFilterMultimapNmax 20 \
	     --alignSJoverhangMin 8 \
	     --alignSJDBoverhangMin 1 \
	     --outFilterMismatchNmax 999 \
	     --outFilterMismatchNoverReadLmax 0.04 \
	     --alignIntronMin 20 \
	     --alignIntronMax 1000000 \
	     --alignMatesGapMax 1000000 \
	     > $STAR_loc/${experiment}_star.log 2>&1

else #PE

	echo "Aligning with STAR (PE) - logs will be in a seperate log file"
	STAR --runThreadN $numbproc \
	     --genomeDir "$genome" \
	     --readFilesIn "$r1" "$r2"\
	     --readFilesCommand zcat \
	     --outFileNamePrefix $STAR_loc \
	     --outSAMtype BAM SortedByCoordinate \
	     --outSAMunmapped Within \
	     --outFilterType BySJout \
	     --outFilterMultimapNmax 20 \
	     --alignSJoverhangMin 8 \
	     --alignSJDBoverhangMin 1 \
	     --outFilterMismatchNmax 999 \
	     --outFilterMismatchNoverReadLmax 0.04 \
	     --alignIntronMin 20 \
	     --alignIntronMax 1000000 \
	     --alignMatesGapMax 1000000 \
	     > $STAR_loc/${experiment}_star.log 2>&1
fi


#Move the STAR results to intended - hard drive - location.
mv $STAR_loc $out
samtools index "$out/STAR_res/Aligned.sortedByCoord.out.bam"


