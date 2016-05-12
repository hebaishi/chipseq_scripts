#!/bin/sh
columnCount() {
  return `head -n1 $1 | sed -e 's|\s\+|\n|g' |wc -l`
}

checkWarn() {
  if [ -z `which $1`  ]; then
    echo Error: $1 not found in your path
    return 0
  else
    return 1
  fi
}

if [ "$#" -ne 2 ]; then
    echo Usage: overlap_peaks.sh file1.bed file2.bed
    echo
    echo overlap_peaks.sh accepts as inputs two BED files that
    echo you wish to overlap/intersect.
    echo The output is a tab-delimited text file without headers.
    echo For each overlap, the script prints the relevant row from
    echo the first BED file then the second BED file, followed by
    echo the coordinates and length of the core overlap, followed
    echo by the coordinates and length of the larger, merged interval.
    exit
fi

checkWarn bedtools
if [ "$?" -lt "1" ]; then  exit; fi

columnCount $1
file_one_cols=$?

bedtools intersect -wa -wb -a $1 -b $2 | awk "OFS = \"\\t\"{
  if ( \$2 < \$(2 + $file_one_cols) )
  {merge_start = \$2; intersect_start = \$(2 + $file_one_cols); } else {
    intersect_start = \$2; merge_start = \$(2 + $file_one_cols); }

  if ( \$3 < \$(3 + $file_one_cols) )
  { intersect_end = \$3; merge_end = \$(3 + $file_one_cols); } else {
    merge_end = \$3; intersect_end = \$(3 + $file_one_cols);}

    print \$0 \"\\t\" intersect_start \"\\t\" intersect_end \"\\t\" (intersect_end - intersect_start \"\\t\" merge_start \"\\t\" merge_end \"\\t\" (merge_end - merge_start) )
}"
