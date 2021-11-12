#!/bin/bash

while getopts 'c:d:g:s:i:w:' c
do
  case $c in
    c) CONTIGS=$OPTARG ;;
    d) DOMAIN=$OPTARG ;;
    g) GCODE=$OPTARG ;;
    s) SPECIES=$OPTARG ;;
    i) ISOLATE=$OPTARG ;;
    w) WORKFLOW=$OPTARG ;;
  esac
done

# Set printed-to-screen script text color (generally don't change):
COLOR='\033[1;36m'
NC='\033[0m'

INPUT_GTO=${ISOLATE}.input.gto
OUTPUT_GTO=${ISOLATE}.output.gto
OUTPUT_GFF=${ISOLATE}.gff
OUTPUT_GBK=${ISOLATE}.gbk

printf "${COLOR}\nCreating rasttk genome object...${NC}\n"
rast-create-genome --scientific-name ${SPECIES} --genetic-code ${GCODE} --domain ${DOMAIN} --contigs ${CONTIGS} > ${INPUT_GTO}

printf "${COLOR}\nAnnotating genome...${NC}\n"
rast-process-genome --workflow ${WORKFLOW} < ${INPUT_GTO} > ${OUTPUT_GTO}

printf "${COLOR}\nExporting GFF file...${NC}\n"
rast-export-genome gff < ${OUTPUT_GTO} > ${OUTPUT_GFF}

printf "${COLOR}\nExporting Genbank file...${NC}\n"
rast-export-genome genbank < ${OUTPUT_GTO} > ${OUTPUT_GBK}

printf "${COLOR}\nDone.${NC}\n"

