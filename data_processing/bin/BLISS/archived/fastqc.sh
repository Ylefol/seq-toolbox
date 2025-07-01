#!/usr/bin/env bash

# Input file and base name
DIR="$1"
FILE="$2"
BASE=$(basename "$FILE" .fastq.gz)

# Output directory
OUTDIR="$3"

# Define the path to the FastQC executable
FASTQC=$(which fastqc)

# Create an error log file
FASTQC_ERR="$OUTDIR/$BASE.fastqc.err"

FILE="${FILE}.fastq.gz"

# Print a message to the log
echo "Running FastQC: $BASE"

# Run FastQC
fastqc "$DIR/$FILE" --outdir "$OUTDIR" 2> "$FASTQC_ERR"

#Fastq screen
fastq_screen --conf /media/yohanl/YOHAN_scTiSA/fastq_screen_conf/fastq_screen.conf "$DIR/$FILE" --aligner Bowtie2 --outdir "$OUTDIR"


# Check if FastQC completed successfully
if [ $? -eq 0 ]; then
    echo "FastQC completed successfully for: $BASE"
else
    echo "FastQC encountered errors. Check the log: $FASTQC_ERR"
fi
echo 'Done'

