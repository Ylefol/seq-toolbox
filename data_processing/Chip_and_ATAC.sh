#!/usr/bin/env bash

################################################################################
# clear
# DEFINING VARIABLES
experiment=$1			# experiment ID found in fastq filename: expID_R1.fastq.gz
genome=$2			# human or mouse
fastqDir=$3			# full/path/to/directory containing the compressed fastq file
numbproc=$4			# set the desired number of threads to be used while mapping to reference genome
################################################################################

# PREPARE DIRECTORY STRUCTURE
datadir=./results && mkdir -p $datadir/$experiment
finalRes=$datadir/$experiment/peaks_bw_bed && mkdir -p $finalRes
out=$datadir/$experiment/outdata && mkdir -p $out
aux=$datadir/$experiment/auxdata && mkdir -p $aux
# create a folder for QC results
QC_out=$datadir/$experiment/QC_results && mkdir -p $QC_out

#May have to be automated
script_dir="/home/yohanl/A_Projects/seq-toolbox/data_processing/"

# Create some parameters - primarily checking file names to see if it will be a single-end or paired-end
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

#Set start time
start=`date +%s`

#set genome references, blacklists, and effectivegenomesize
################################################################################
read blacklisted_region effective_genome_size ref_genome < <("$script_dir"bin/set_refgen.sh $genome)

#FASTQC and FAST_SCREEN
################################################################################

#"$script_dir"bin/fastqc_fastscreen.sh $r1 $QC_out
#if [ $numb_of_files == 2 ]; then
#	"$script_dir"bin/fastqc_fastscreen.sh $r2 $QC_out
#fi


#TRIMMOMATIC SE AND PE
################################################################################
#trim_file="All_adapters.fa"
#trim_file="$script_dir"trimmomatic_files/"$trim_file"

#"$script_dir"bin/trimmomatic.sh $r1 $r2 $out $experiment $trim_file $QC_out

#trim_1=$out/$experiment"_trimmed.fastq.gz"
#trim_2="NONE" #Temporary name allocation
#if [ $numb_of_files == 2 ]; then
#	trim_1=$out/$experiment"_trimmed_paired_r1.fastq.gz"
#	trim_2=$out/$experiment"_trimmed_paired_r2.fastq.gz"
#fi


#Mapping with Bowtie2 and conversion to BAM
################################################################################
trim_1="/media/yohanl/Expansion/seq-processing/Lisa_5hmu/results_T2T"/$experiment/"outdata"/$experiment"_trimmed_paired_r1.fastq.gz"
trim_2="/media/yohanl/Expansion/seq-processing/Lisa_5hmu/results_T2T"/$experiment/"outdata"/$experiment"_trimmed_paired_r2.fastq.gz"
"$script_dir"bin/bowtie2.sh $trim_1 $trim_2 $experiment $out $ref_genome $numbproc


#BLACKLISTING
################################################################################
echo "################################################################################"

if [[ "$your_variable" != "None" ]]; then
    echo "blacklist found: $blacklisted_region" 
    bedtools intersect -v -abam $out/$experiment"_sorted.bam" -b "$blacklisted_region" | samtools sort -@ "$numbproc" - -o $out/$experiment"_sorted.bam"
else
echo "no blacklist found"
fi

echo "################################################################################"


#Indexing and extra filtering
################################################################################

"$script_dir"bin/processing_and_filtering.sh $out $aux $experiment $numbproc


#Run bamQC on the bam files before and after filters
################################################################################
echo "################################################################################"
echo "Running BAMQC"
echo "################################################################################"
echo "BAMQC for sorted results"
qualimap bamqc -bam $out/$experiment"_sorted.bam" -outdir $QC_out -outfile $experiment"_bamQC_sorted.pdf"
echo "################################################################################"
echo "BAMQC for dedup results"
qualimap bamqc -bam $out/$experiment"_dedup.bam" -outdir $QC_out -outfile $experiment"_bamQC_dedup.pdf"


#Create bigwig, bed files, and MACS3 peak calling
################################################################################
echo "################################################################################"
echo "Creating basic results"
echo "################################################################################"
"$script_dir"bin/produce_basic_results.sh $out $experiment $finalRes $effective_genome_size $numbproc


#Remove auxilliary folder and files
rm -rf $aux

end=`date +%s`

runtime=$((end-start))

echo "################################################################################"
echo "Pipeline completed in $runtime seconds"
echo "################################################################################"






