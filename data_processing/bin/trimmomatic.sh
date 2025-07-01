#!/usr/bin/env bash

r1=$1
r2=$2
output_loc=$3
output_name=$4
trim_file=$5
QC_dir=$6

trim_param="ILLUMINACLIP:$trim_file:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36"
if [ "$r2" = "NONE" ]; then
	
	echo "################################################################################"
	echo "Running Trimmomatic (SE) for: $output_name"
	echo "################################################################################"
	trimmomatic SE -phred33 $r1 $output_loc/$output_name"_trimmed.fastq.gz" $trim_param
	echo "################################################################################"
	echo "Running fastqc for trimmed file"
	echo "################################################################################"
	fastqc $output_loc/$output_name"_trimmed.fastq.gz" --outdir $QC_dir
 
else

	echo "################################################################################"
	echo "Running Trimmomatic (PE) for: $output_name"
	echo "################################################################################"
	trimmomatic PE -phred33 $r1 $r2 $output_loc/$output_name"_trimmed_paired_r1.fastq.gz"  $output_loc/$output_name"_trimmed_unpaired_r1.fastq.gz" $output_loc/$output_name"_trimmed_paired_r2.fastq.gz" $output_loc/$output_name"_trimmed_unpaired_r2.fastq.gz" $trim_param
	echo "################################################################################"
	echo "Running fastqc for trimmed files"
	echo "################################################################################"
	fastqc $output_loc/$output_name"_trimmed_paired_r1.fastq.gz" --outdir $QC_dir
	fastqc $output_loc/$output_name"_trimmed_paired_r2.fastq.gz" --outdir $QC_dir

fi


