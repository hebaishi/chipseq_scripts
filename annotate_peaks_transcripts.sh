#!/bin/bash
versionSplit() {
  return `echo $1 | sed -e 's|\.|\n|g' | head -n$2 | tail -n1`
}

check_version_greater_equal() {
  version_number_current=`$2 2>&1 | grep -oP "\d+\.\d+\.\d+"`
  version_number_expected=`echo $3 | grep -oP "\d+\.\d+\.\d+"`
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
    exit
  fi
}

checkSorted() {
  if [ `sort -c -k1,1 -k2,2n $1 2>&1 | wc -l` -gt "0" ]; then
    echo File $1 is NOT sorted!
    echo Please re-run with a sorted BED file!
    exit
  fi
}

if [ "$#" -ne 4 ]; then
  echo
  echo Overview:
  echo =========
  echo annotate_peaks_transcripts.sh accepts as inputs a BED file of your favourite
  echo peaks/intervals, and a standard Ensembl GTF file. The output is a
  echo tab-delimited text file with headers.
  echo
  echo This can be run in two modes: genes within a specified
  echo distance of peaks, or closest gene\(s\) to peaks.
  echo
  echo dist mode:
  echo ==========
  echo Usage: annotate_peaks_transcripts.sh dist file.bed file.gtf distance_in_bp
  echo
  echo In \"dist\" mode, the script calculates peak centers and TSS-peak center
  echo distances. As the script uses the entirety of a peak to calculate distance
  echo to TSSs, please note that there may be cases where the TSS-peak center
  echo distance is greater than the user-specified distance.
  echo
  echo clos mode:
  echo ==========
  echo Usage: annotate_peaks_transcripts.sh clos file.bed file.gtf N_closest
  echo
  echo In \"clos\" mode, the script calculates the smallest distance between the
  echo peak \(unstranded\) and TSS. In the case of an overlap, note that the script
  echo returns a zero value for this distance. Note also that the script returns
  echo \".\" or \"-1\" in genomic feature fields when a peak does not coincide with
  echo any transcripts in a chromosome/scaffold.
  exit
fi

checkWarn gtf2tab
checkWarn bedtools

check_version_greater_equal 'gtf2tab' 'gtf2tab' v0.1.20
check_version_greater_equal 'bedtools' 'bedtools --version' v2.23.0

# Create temporary file
temp_annotation_file=`mktemp`

if [ "$1" = "dist" ]; then
  # Generate temporary annotation file
  gtf2tab -C 1 -t transcript -f 1,4,5,7 -a gene_id,transcript_id,gene_name,gene_biotype $3 | awk "OFS=\"\\t\" { if ( \$4 == \"+\") {start_pos =\$2 } else {start_pos = \$3}; \$2 = start_pos - $4; if (\$2 < 0) {\$2 = 0;} \$3 = start_pos +  $4; \$9 = start_pos; print  }" | tail -n+2 > $temp_annotation_file
  # Print header
  echo -e "chr\tTSS_minus_$3\tTSS_plus_$3\tstrand\tgene_id\ttranscript_id\tgene_name\tgene_biotype\tTSS\tpeak_center\tTSS_to_peak_center\tbed_chr\tbed_start\tbed_end\tbed_other_columns"
  # Perform the intersection and calculate peak center + distance to TSS
  bedtools intersect -wa -wb -a $temp_annotation_file -b $2 | awk 'OFS = "\t" { peak_center = int( (($12 - $11)/ 2) + $11  ); $9 = $9"\t" peak_center "\t" (peak_center - $9); print;  }'
elif [ "$1" = "clos" ]
then
  checkSorted $2
  # Generate temporary annotation file
  gtf2tab -C 1 -t transcript -f 1,4,5,7 -a gene_id,transcript_id,gene_name,gene_biotype $3 | awk "OFS=\"\\t\" { if ( \$4 == \"+\") {start_pos =\$2 } else {start_pos = \$3}; \$2 = start_pos; \$3 = start_pos; print }" | tail -n+2 | sort -k1,1 -k2,2n > $temp_annotation_file
  bedtools closest -t all -k $4 -D b -b $temp_annotation_file -a $2
fi

rm $temp_annotation_file
