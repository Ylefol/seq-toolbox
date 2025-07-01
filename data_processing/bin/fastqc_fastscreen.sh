#!/usr/bin/env bash


file=$1
QC_loc=$2


# Create an error log
echo "################################################################################"
echo "Running fastqc for: $file"
echo "################################################################################"
fastqc $file --outdir "$QC_loc"

#Fastq screen
echo "################################################################################"
echo "Running fastq_screen for: $file"
echo "################################################################################"
fastq_screen --conf /media/yohanl/Expansion/reference_genomes/bowtie2/fastq_screen.conf $file --aligner Bowtie2 --outdir "$QC_loc"

