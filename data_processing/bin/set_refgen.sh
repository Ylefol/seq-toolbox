#!/usr/bin/env bash

genome=$1
aligner=$2

#set genome references, blacklists, and effectivegenomesize
################################################################################
if [[ $genome == GRCh38 ]]; then
blacklisted_region="/home/yohanl/A_Projects/seq-toolbox/data_processing/blacklist_files/hg38-blacklist.v2.bed.gz"
effective_genome_size="2913022398"
ref_genome="/media/yohanl/Expansion/reference_genomes/"$aligner"/GRCh38/GRCh38"
fi
if [[ $genome == GRCm38 ]]; then
blacklisted_region="/home/yohanl/A_Projects/seq-toolbox/data_processing/blacklist_files/mm10-blacklist.v2.bed.gz"
effective_genome_size="2652783500"
ref_genome="/media/yohanl/Expansion/reference_genomes/bowtie2/GRCm38/GRCm38"
fi
if [[ $genome == WBcel235 ]]; then
blacklisted_region="/home/yohanl/A_Projects/seq-toolbox/data_processing/blacklist_files/ce11-blacklist.v2.bed.gz"
effective_genome_size="100272607" #Calculated using faCount on ce11 fasta file - expected to be a bit over but good enough
ref_genome="/media/yohanl/Expansion/reference_genomes/bowtie2/WBcel235/WBcel235"
fi
if [[ $genome == BDGP6 ]]; then
blacklisted_region="/home/yohanl/A_Projects/seq-toolbox/data_processing/blacklist_files/dm6-blacklist.v2.bed.gz"
effective_genome_size="142573017"
ref_genome="/media/yohanl/Expansion/reference_genomes/bowtie2/BDGP6/BDGP6"
fi
if [[ $genome == T2TCHM13 ]]; then
blacklisted_region="/home/yohanl/A_Projects/seq-toolbox/data_processing/blacklist_files/T2TCHM13-blacklist.bed.gz"
effective_genome_size="3117292070"
ref_genome="/media/yohanl/Expansion/reference_genomes/bowtie2/T2TCHM13/T2TCHM13v2.0"
fi

echo $blacklisted_region $effective_genome_size $ref_genome
