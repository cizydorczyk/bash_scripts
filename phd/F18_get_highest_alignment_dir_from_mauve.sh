#!/bin/bash

# $1 = parent directory
# $2 = output directory
# $3 = isolate

shopt -s extglob
DIR=$1 HIGHEST=''
for A in "$DIR"/alignment+([[:digit:]]); do
    B=${A#*/alignment}
    [[ -z $HIGHEST || B -gt HIGHEST ]] && HIGHEST=$B
done
[[ -n $HIGHEST ]] && echo "Highest: $DIR/alignment${HIGHEST}"

cp $DIR"/alignment"${HIGHEST}"/"$3"_filtered_assembly.fasta" $2

mv $2"/"$3"_filtered_assembly.fasta" $2"/"$3"_filtered_ordered_assembly.fasta"
