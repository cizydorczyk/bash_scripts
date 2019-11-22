#!/bin/bash

# Variables:
# 1: snippy-core output directory (where all output will be generated)
# 2: snippy-core output prefix

# Set printed-to-screen script text color (generally don't change):
COLOR='\033[1;36m'
NC='\033[0m'

cd $1

printf "${COLOR}Removing reference genome from whole genome alignment...${NC}\n"

python ~/python_scripts/phd/snippy-remove-ref.py --aln $2".full.aln" --out "no_ref."$2".full.aln"

printf "${COLOR}Cleaning no ref. whole genome alignment...${NC}\n"

snippy-clean_full_aln "no_ref."$2".full.aln" > "clean.no_ref."$2".full.aln"

printf "${COLOR}Removing uncleaned no ref. whole genome alignment...${NC}\n"

rm "no_ref."$2".full.aln"

printf "${COLOR}Removing ref. from SNP alignment...${NC}\n"

python ~/python_scripts/phd/snippy-remove-ref.py --aln $2".aln" --out "no_ref.intermediate."$2".aln"

printf "${COLOR}Removing any identical positions from no ref. SNP alignment...${NC}\n"

snp-sites -c -o "snps-only.no_ref."$2".aln" "no_ref.intermediate."$2".aln"

printf "${COLOR}Removing intermediate SNP alignment with possible identical sites due to removing ref...${NC}\n"

rm "no_ref.intermediate."$2".aln"

printf "${COLOR}Generating matrix of SNP distances before removing recombination, etc. ...${NC}\n"

snp-dists "snps-only.no_ref."$2".aln" > "snps-only.no_ref."$2".aln.matrix.txt"

printf "${COLOR}Done.${NC}\n"
