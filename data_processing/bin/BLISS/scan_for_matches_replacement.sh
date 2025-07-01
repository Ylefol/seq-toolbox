#!/usr/bin/env bash

input="$1"  #input file - full path
output="$2"  #output file - full path
barcode="$3"  # sample barcode
max_mismatches=$4  # number of allowed mismatches


cat $input | parallel --block 100M -k --pipe -L 2 '
    awk -v barcode="'$barcode'" -v max_mismatches='$max_mismatches' "
    function hamming_distance(s1, s2) {
        dist = 0
        for (i = 1; i <= length(s1); i++) {
            if (substr(s1, i, 1) != substr(s2, i, 1)) {
                dist++
                if (dist > max_mismatches) return dist
            }
        }
        return dist
    }

    NR % 2 == 1 {header = \$0; next}  # Store header

    {
        sub_seq = substr(\$0, 9, 8)  # Extract nucleotides 9-16
        if (hamming_distance(sub_seq, barcode) <= max_mismatches) {
            seq_length = length(\$0)
            formatted_seq = \" \" substr(\$0,1,8) \" \" substr(\$0,9,8) \" \" substr(\$0,17) \" \" \" \"  # Add spaces
            print header \":\" \"[1,\" seq_length \"]\"
            print formatted_seq
        }
    }"
' > $output
