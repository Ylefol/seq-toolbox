#!/usr/bin/env bash

################################################################################
# DEFINING VARIABLES
experiment=$1			# experiment ID found in fastq filename: expID_R1.fastq.gz
genome=$2			# human or mouse
fastqDir=$3			# full/path/to/directory containing the compressed fastq file
numbproc=$4			# set the desired number of threads to be used while mapping to 

# PREPARE DIRECTORY STRUCTURE
datadir=./results && mkdir -p $datadir/$experiment
out=$datadir/$experiment/outdata && mkdir -p $out
# create a folder for QC results
QC_out=$datadir/$experiment/QC_results && mkdir -p $QC_out

#May have to be automated
script_dir="/home/yohanl/A_Projects/seq-toolbox/data_processing/"

# Create some parameters - primarily checking file names to see if it will be a single-end or paired-end
################################################################################
find $fastqDir -maxdepth 1 -type f -iname "${experiment}*.fastq.gz" | sort > filelist_"$experiment"

numb_of_files=`cat filelist_"$experiment" | wc -l`
r1=`cat filelist_"$experiment" | head -n1`
r2="NONE" #Create and empty R2, will be overwritten if a real one exists
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

#FASTQC and FAST_SCREEN
################################################################################

#"$script_dir"bin/fastqc_fastscreen.sh $r1 $QC_out
#if [ $numb_of_files == 2 ]; then
#	"$script_dir"bin/fastqc_fastscreen.sh $r2 $QC_out
#fi

#TRIMMOMATIC SE AND PE
################################################################################
trim_file="All_adapters.fa"
trim_file="$script_dir"trimmomatic_files/"$trim_file"

#"$script_dir"bin/trimmomatic.sh $r1 $r2 $out $experiment $trim_file $QC_out


#FASTQC OF TRIMMED FILES
################################################################################
trim_1=$out/$experiment"_trimmed.fastq.gz"
trim_2="NONE" #Temporary name allocation
if [ $numb_of_files == 2 ]; then
	trim_1=$out/$experiment"_trimmed_paired_r1.fastq.gz"
	trim_2=$out/$experiment"_trimmed_paired_r2.fastq.gz"
fi


# Run STAR and capture logs
################################################################################

"$script_dir"bin/STAR.sh $trim_1 $trim_2 $out $experiment $genome $numbproc $script_dir


# Run FeatureCounts -- start file here
################################################################################
"$script_dir"bin/featureCounts.sh $out $experiment $genome $numbproc 

#End feature counts file here

end=`date +%s`

runtime=$((end-start))

echo "################################################################################"
echo "Pipeline completed in $runtime seconds"
echo "################################################################################"






