#! /usr/bin/env bash

input=$1
genome=$2
output=$3

if [ $genome == human ]; then
    cat $input |
    #This performs some sort of adjustment, not sure if required, removing.
    #awk -F'\t' '{ if ( $4 == "-" ) {$2=$2-1;$3=$3-1;print $0} else {print $0} }'| #correct for the bedtools/bamfile mismatch on negative strands end location 
    #tr ' ' '\t' |
    grep -v "_" | sed -e 's/chrX/chr23/g' | sed -e 's/chrY/chr24/g' | awk '{if($4=="+"){printf "%s\t%d\t%d\n", $1, $2, $3}else{printf "%s\t%d\t%d\n", $1, $3, $3+1}}' | sort -k1,1 -k2,2g -k3,3g | LC_ALL=C uniq -c | awk '{OFS="\t";print $2,$3,$4,$1}' > $output
fi
if [ $genome == mouse ]; then
    cat $input | grep -v "_" | sed -e 's/chrX/chr21/g' | sed -e 's/chrY/chr22/g' | awk '{if($4=="+"){printf "%s\t%d\t%d\n", $1, $2, $3}else{printf "%s\t%d\t%d\n", $1, $3, $3+1}}' | sort -k1,1 -k2,2g -k3,3g | LC_ALL=C uniq -c | awk '{OFS="\t";print $2,$3,$4,$1}' > $output
fi




