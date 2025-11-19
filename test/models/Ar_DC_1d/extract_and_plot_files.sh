#!/bin/bash

glob_pattern=${1:-"plt?????"}
varnum=${2:-"2"}
#echo "$glob_pattern"
for file in $glob_pattern; do
  #echo "Processing file: $file"
  #assuming amrex's fextract is in
  #your home directory
  ~/fextract "$file"
done

gnuplot -persist -c plot_all.gp "$glob_pattern".slice $varnum
