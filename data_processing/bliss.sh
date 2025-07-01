#!/usr/bin/env bash

################################################################################
# clear
# DEFINING VARIABLES
experiment=$1		# experiment ID found in fastq filename: expID_R1.fastq.gz
genome=$2			# human or mouse
barcode=$3			# is the UMI+barcode pattern file used in the linker
quality=$4			# minimum mapping quality desired
fastqDir=$5			# full/path/to/directory containing the compressed fastq file
numbproc=$6			# set the desired number of threads to be used while mapping to reference genome
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

"$script_dir"bin/fastqc_fastscreen.sh $r1 $QC_out
if [ $numb_of_files == 2 ]; then
	"$script_dir"bin/fastqc_fastscreen.sh $r2 $QC_out
fi


# Mapping and formatting output
################################################################################

# Create a 'oneline' format of the fastq file(s). This allows us to reconstitute a fastq file after we have used scan_for_matches and before those results are mapped
"$script_dir"bin/BLISS/prepare_files.sh  $r1 $in $numb_of_files $r2

echo "Filtering reads based on patterns ..."
#cat $in/r1.fa | parallel --block 100M -k --pipe -L 2 "scan_for_matches $patfile - " > $out/filtered.r1.fa
#A replacement for scan_for_matches - parameters are input, output, barcode, then number of allowed mismatches
"$script_dir"bin/BLISS/scan_for_matches_replacement.sh $in/r1.fa $out/filtered.r1.fa $barcode 1
echo "Done"
# This reconstitutes a fastq file for mapping. It requires the 'one line' format of the files as well as the results of scan_for_matches
"$script_dir"bin/BLISS/prepare_for_mapping.sh $numb_of_files $out $aux $in

# Aligns to the reference genome provided using bwa
"$script_dir"bin/BLISS/mapping.sh $numb_of_files $numbproc $refgen $aux $out $experiment $quality

# Here we join and format results. We join the mapped results with the scan_for_matches results. These are then formatted
"$script_dir"bin/BLISS/umi_joining.sh $numb_of_files $out $experiment $aux $quality

# Further formatting of the outputted bed file
cat $out/_q"$quality".bed | cut -f-5 |LC_ALL=C uniq -c | awk '{print $2,$3,$4,$5,$6,$1}' | tr " " "," > $out/pre_umi_filtering.csv

#UMI filtering
################################################################################

# Group consecutive UMIs as well as UMIs which are spatially close to each other (space gap of 30 mm gap of 2)
python3 "$script_dir"bin/BLISS/umi_filtering.py $out/pre_umi_filtering.csv $out/q$quality"_aux"

#Formatting - remove commas, apostrophe, brackets (left and right) and replace spaces with tabs (\t)
cat $out/q"$quality"_aux | tr -d "," | tr -d "'" | tr -d "[" | tr -d "]" | tr " " "\t" > $out/q"$quality"_chr-loc-strand-umi-pcr

# Format the results. This formatting renames chromosomes, for example in human chrx and chry become chr 23 and chr24 respectively.
# This changes from one organism to another therefore the genome is also included in the command. 
"$script_dir"bin/BLISS/umi_filter.sh $out/q"$quality"_chr-loc-strand-umi-pcr $genome $out/q"$quality"_chr-loc-countDifferentUMI.bed

#Create summary statistics.
################################################################################

echo "Alignment statistics:" >> $out/summary.txt
samtools flagstat --threads $numbproc $out/*.all.bam >> $out/summary.txt

echo "Number of reads on plus strand and on minus strand:" >> $out/summary.txt
cat $out/q"$quality"_chr-loc-strand-umi-pcr | grep -v "_" | cut -f4 | sort | uniq -c >> $out/summary.txt
echo "Number of DSB locations:" >> $out/summary.txt
cat $out/q"$quality"_chr-loc-countDifferentUMI.bed | grep -v "_" | wc -l >> $out/summary.txt
echo "Number of UMIs:" >> $out/summary.txt
cat $out/q"$quality"_chr-loc-strand-umi-pcr | grep -v "_" | wc -l >> $out/summary.txt

name=`echo $patfile|rev|cut -d'/' -f1|rev`
cat $out/q"$quality"_chr-loc-strand-umi-pcr |
    awk -F'\t' '{ if ( $4 == "-" ) {$2=$2-1;$3=$3-1;print $0} else {print $0} }'| #correct for the bedtools/bamfile mismatch on negative strands end location 
    tr ' ' '\t' > $out/chr-loc-strand-umi-pcr.tsv
#mv $out/q"$quality"_chr-loc-countDifferentUMI.bed $out/chr-loc-countDifferentUMI.bed

#Clean and format directory
################################################################################
#Delete unneeded folders and files
rm -r "$in"*
rm -r "$aux"*
rm $out/filtered.r1.fa
rm $out/_q60.bed
rm $out/q60_aux

#Move results from outdata folder to results folder (outdata is contained within results but it no longer necessary at the end of the pipeline run)
mv "$out"/* $datadir/$experiment/
rm -r "$out"* 

#Bit of extra sorting
extra_out=$datadir/$experiment/extra_files && mkdir -p $extra_out

mv $datadir/$experiment/pre_umi_filtering.csv $extra_out/pre_umi_filtering.csv
mv $datadir/$experiment/q"$quality"_chr-loc-strand-umi-pcr $extra_out/q"$quality"_chr-loc-strand-umi-pcr
mv $datadir/$experiment/"$experiment".all.bam $extra_out/"$experiment".all.bam
mv $datadir/$experiment/"$experiment".all.bam.bai $extra_out/"$experiment".all.bam.bai
mv $datadir/$experiment/"$experiment".q"$quality".bam $extra_out/"$experiment".q"$quality".bam
mv $datadir/$experiment/"$experiment".q"$quality".bam.bai $extra_out/"$experiment".q"$quality".bam.bai



