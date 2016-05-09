#!/bin/bash

versionSplit() {
  return `echo $1 | sed -e 's|\.|\n|g' | head -n$2 | tail -n1`
}

check_version_greater_equal() {
  version_number_current=`$1 2>&1 | grep -oP "\d+\.\d+\.\d+"`
  version_number_expected=`echo $2 | grep -oP "\d+\.\d+\.\d+"`
  for idx in $(seq 1 3);
    do versionSplit $version_number_current $idx
    version_current=$?
    versionSplit $version_number_expected $idx
    version_expected=$?
    if [ "$version_expected" -gt "$version_current" ]
    then
      echo Error: Your version of $1 [ $version_number_current ] is too old!
      echo Please update $1 to continue
      return 0
    elif [ "$version_expected" -lt "$version_current" ]
    then
      return 1
    fi
  done
  return 1
}

checkWarn() {
  if [ -z `which $1`  ]; then
    echo Error: $1 not found in your path
    return 0
  else
    return 1
  fi
}

if [ "$#" -ne 3 ]; then
    echo Usage: overlap_peaks_genes.sh file.bed file.gtf distance_in_bp
    exit
fi

checkWarn gtf2tab
if [ "$?" -lt "1" ]; then  exit; fi

checkWarn bedtools
if [ "$?" -lt "1" ]; then  exit; fi

check_version_greater_equal 'gtf2tab' v0.1.19
if [ "$?" -lt "1" ]; then  exit; fi

# Create temporary file
temp_annotation_file=`mktemp`
# Generate temporary annotation file
gtf2tab -C 1 -t transcript -f 1,4,5,7 -a gene_id,transcript_id,gene_name,gene_biotype $2 | awk "OFS=\"\\t\" { if ( \$4 == \"+\") {start_pos =\$2 } else {start_pos = \$3}; \$2 = start_pos - $3; if (\$2 < 0) {\$2 = 0;} \$3 = start_pos +  $3; \$9 = start_pos; print  }" | tail -n+2 > $temp_annotation_file
# Print header
echo -e "chr\tTSS_minus_$3\tTSS_plus_$3\tstrand\tgene_id\ttranscript_id\tgene_name\tgene_biotype\tTSS\tpeak_center\tTSS_to_peak_center\tbed_chr\tbed_start\tbed_end\tbed_other_columns"
# Perform the intersection and calculate peak center + distance to TSS
bedtools intersect -wa -wb -a $temp_annotation_file -b $1 | awk 'OFS = "\t" { peak_center = int( (($12 - $11)/ 2) + $11  ); $9 = $9"\t" peak_center "\t" (peak_center - $9); print;  }'
# Delete temporary file
rm $temp_annotation_file
