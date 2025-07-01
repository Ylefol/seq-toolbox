#!/bin/usr/env bash

#Set start time
start=`date +%s`

#!/usr/bin/env bash
# save as summarize_barcodes.sh and run with: bash summarize_barcodes.sh /path/to/dir

DIR="$1"
OUTFILE="barcode_counts.txt"

# write header
echo -e "filename\tcount\tbarcode" > "$OUTFILE"

# iterate all .fastq.gz files
for f in "$DIR"/*.fastq.gz; do
  # skip if no files match
  [[ -e "$f" ]] || continue
  echo "Checking $f"
  # compute most‐common 8‐mer (pos 8–15) and its count
  read count barcode < <(
    zcat "$f" \
      | awk 'NR%4==2 { c=substr($0,8,8); counts[c]++ }
             END { for (k in counts) print counts[k], k }' \
      | sort -k1,1nr \
      | head -n1
  )

  # append filename (basename), count and barcode
  echo -e "$(basename "$f")\t${count:-0}\t${barcode:-NA}" >> "$OUTFILE"
done


echo "Took $runtime seconds"



