#!/bin/bash
checkWarn() {
  if [ -z `which $1`  ]; then
    echo Error: $1 not found in your path
    exit
  fi
}

checkWarn bedtools

if [ "$#" -ne 4 ]; then
    echo Usage: weeder_sequence_extract.sh [gtf_file] [chromosome_sizes_file] [genome_fasta_file] [distance_in_bp]
    exit
fi

gtf2tab -t transcript -f 1,4,5,7 $1 |tail -n+2 | awk 'OFS="\t"{strand=$4; $4="peak"NR; $5=0; $6=strand; print;}' | bedtools flank -l $4 -r 0 -s -i /dev/stdin -g $2 | bedtools getfasta -fi $3 -bed /dev/stdin -fo /dev/stdout
