#!/usr/bin/env bash
numb_of_files=$1
out=$2
sample=$3
aux=$4
quality=$5

echo 'Selecting unique UMIs'
if [[ $numb_of_files == 1 ]]; then # IF SE READS
    # Convert bam to bed - first awk select lines where strand is + then customize the output for that strand and save is as forward
    bedtools bamtobed -i ${out}/${sample}.q${quality}.bam | awk '$6 == "+"' | awk '{OFS="\t";print $4,$1,$2,"+"}' > $aux/forward & pid1=$! # if + strand DSB location is the second field
    bedtools bamtobed -i ${out}/${sample}.q${quality}.bam | awk '$6 == "-"' | awk '{OFS="\t";print $4,$1,$3,"-"}' > $aux/reverse & pid2=$! # if - strand DSB location is the third field
    
    # These commands ensure that both of the above files are done before proceeding
    wait $pid1
    wait $pid2
    
    # combines forward and reverse file into a single file
    cat $aux/forward $aux/reverse | LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -k1,1 > $aux/id.chr.loc.strand & pid1=$!
    # Takes the scan_for_matches file, removes '>' and keeps only the first seven fields, pastes everything together and saves in aux
    cat $out/filtered.r1.fa | tr -d ">" | cut -d':' -f-7 | paste - - > $aux/id.umi.barcode.genomic_aux & pid2=$!
    
    #Ensure that both of the above processes are complete
    wait $pid1
    wait $pid2
    
    # clean columns, if there are five fields keep only the relevant ones
    cat $aux/id.umi.barcode.genomic_aux | awk '{if (NF == 4) print; else if (NF == 5) print $1,$3,$4,$5}' | LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -k1,1 > $aux/id.umi.barcode.genomic
    
    # Join strand data with genomic information and format the output
    LC_ALL=C join $aux/id.chr.loc.strand $aux/id.umi.barcode.genomic | cut -d' ' -f-5 |LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -t' ' -k2,2 -k3,3n -k4,4 -k5,5 |
    awk '{print $2,$3,$3+1,$4,$5,$1}' | tr " " "\t" > $out/_q$quality.bed 
fi

#IMPORTANT -- in the case of PE, the joining still uses R1 (and not R2) of scan_to_matches.
# In fact scan_for_matches is only run on one of the two files (if PE). Since PE is simply a sequencing of both ends of the strand it is assumed
# that the results of scan_for_matches will be applicable to both PE files.
if [[ $numb_of_files == 2 ]]; then # IF PE READS
    bedtools bamtobed -i ${out}/${sample}.q${quality}.bam | awk '$6 == "+"' | awk '{OFS="\t";print $4,$1,$2,"+"}' | grep -v '\\2' | sed 's/\/1//' > $aux/forward & pid1=$!
    bedtools bamtobed -i ${out}/${sample}.q${quality}.bam | awk '$6 == "-"' | awk '{OFS="\t";print $4,$1,$3,"-"}' | grep -v '\\2' | sed 's/\/1//' > $aux/reverse & pid2=$!
    wait $pid1
    wait $pid2
    cat $aux/forward $aux/reverse |LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -k1,1 > $aux/id.chr.loc.strand & pid1=$!
    cat $out/filtered.r1.fa | tr -d ">" | cut -d':' -f-7 | paste - - > $aux/id.umi.barcode.genomic_aux & pid2=$!
    wait $pid1
    wait $pid2
    cat $aux/id.umi.barcode.genomic_aux | awk '{if (NF == 4) print; else if (NF == 5) print $1,$3,$4,$5}' | LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -k1,1 > $aux/id.umi.barcode.genomic
    LC_ALL=C join $aux/id.chr.loc.strand $aux/id.umi.barcode.genomic | cut -d' ' -f-5 | LC_ALL=C sort --parallel=8 --temporary-directory=$HOME/tmp -t' ' -k2,2 -k3,3n -k4,4 -k5,5 | 
    awk '{print $2,$3,$3+1,$4,$5,$1}' | tr " " "\t" > $out/_q$quality.bed
fi

echo 'Done'
