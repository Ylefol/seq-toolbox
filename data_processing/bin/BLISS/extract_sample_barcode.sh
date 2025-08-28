#!/bin/usr/env bash

#Set start time
start=`date +%s`

DIR="$1"
OUTFILE=$DIR/"barcode_counts.txt"

#Set start time
start=`date +%s`

echo "Running sample barcode finder, results will be stored in the provided data directory"

# write header
echo -e "filename\tcount\tbarcode" > "$OUTFILE"

# iterate all .fastq.gz files
for f in "$DIR"/*.fastq.gz; do
  # skip if no files match
  [[ -e "$f" ]] || continue

  base=$(basename "$f")

  # if it ends with _2.fastq.gz, skip it
  if [[ "$base" == *_2.fastq.gz ]]; then
    echo "Skipping $base" >&2
    continue
  fi
  echo "Processing $base"
  # otherwise (ends in _1 or anything else), process it
  read count barcode < <(
    zcat "$f" \
      | awk 'NR%4==2 { c=substr($0,8,8); counts[c]++ }
             END { for (k in counts) print counts[k], k }' \
      | sort -k1,1nr \
      | head -n1
  )

  # append filename, count and barcode (defaulting to 0/NA if none)
  echo -e "${base}\t${count:-0}\t${barcode:-NA}" >> "$OUTFILE"
done

end=`date +%s`

runtime=$((end-start))

echo "Took $runtime seconds"



