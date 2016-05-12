#!/bin/sh
checkWarn() {
  if [ -z `which $1`  ]; then
    echo Error: $1 not found in your path
    exit
  fi
}

checkWarn bedtools
checkWarn bigWigToBedGraph

if [ "$#" -ne 2 ]; then
    echo Usage: bigwig2bed.sh file.bw threshold
    echo The second parameter is the threshold:
    echo everything below the threshold is removed
    exit
fi

bigWigToBedGraph $1 /dev/stdout | awk "{if (\$4>$2) {print;}}" | bedtools sort | bedtools merge
