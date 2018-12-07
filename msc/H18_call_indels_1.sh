#!/bin/bash

########### REMEMBER TO SET REF GENOME LENGTH - 150 bp BELOW!!! #################

# Set tool paths (generally don't change):
SAMTOOLS=/home/conrad/Software/samtools-1.3.1/samtools
PICARD=/home/conrad/Software/picard.jar
GATK=/home/conrad/Software/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
BCFTOOLS=/home/conrad/Software/bcftools-1.3.1/bcftools

# Set project folder:
PF=/home/conrad/Data/primary_project_3/reference_alignments/H18_lesb58_indels

# Set directory with sorted, indexed bam files with read groups added (must do this prior to running this script):
BAMINPUT=/home/conrad/Data/primary_project_3/reference_alignments/H18_lesb58_indels/sorted_indexed_bam_files_with_rg

# Set printed-to-screen script text color (generally don't change):
COLOR='\033[1;36m'
NC='\033[0m'

# Input variables:

# 1 = isolate list (text file, one isolate number per line)

cd $PF/reference

printf "${COLOR}Building reference dictionary...${NC}\n"
java -jar $PICARD CreateSequenceDictionary REFERENCE=reference.fasta OUTPUT="reference.dict"
cd ..

mkdir marked_duplicate_bam_files

cd marked_duplicate_bam_files

mkdir picard_markduplicate_metrics

printf "${COLOR}Marking duplicates...${NC}\n"
for i in $(cat $1); do java -Xmx30g -jar $PICARD MarkDuplicates I=$BAMINPUT/$i"_aln_rg.bam" O=$i"_aln_rg_markeddup.bam" M=picard_markduplicate_metrics/$i"_marked_duplicate_metrics.txt"; done

printf "${COLOR}Indexing marked-duplicate bam files...${NC}\n"
for i in $(cat $1); do $SAMTOOLS index $i"_aln_rg_markeddup.bam"; done

cd ..

mkdir raw_vcf_files

cd raw_vcf_files

printf "${COLOR}Running HaplotypeCaller on indexed marked-duplicate bam files...${NC}\n"
for i in $(cat $1); do java -Xmx30g -jar $GATK -T HaplotypeCaller -R $PF/reference/reference.fasta -I $PF/marked_duplicate_bam_files/$i"_aln_rg_markeddup.bam" -o $i"_raw.vcf" -ploidy 1 -nct 8 -stand_call_conf 20; done

cd ..

mkdir raw_indels

cd raw_indels

printf "${COLOR}Selecting indels only...${NC}\n"
for i in $(cat $1); do java -Xmx30g -jar $GATK -T SelectVariants -R $PF/reference/reference.fasta -nt 8 -V $PF/raw_vcf_files/$i"_raw.vcf" -selectType INDEL -o $i"_raw_indels.vcf"; done

cd ..

mkdir filtered_indel_vcf_files

cd filtered_indel_vcf_files

printf "${COLOR}Filtering indels based on: QUAL>=30, DP>=20, POS>150, POS<6601607...${NC}\n"
################ SET REF LEN - 150 bp IN LINE BELOW WHERE POS<(some number) #######################
for i in $(cat $1); do $BCFTOOLS filter -i 'QUAL>=30 & DP>=20 & POS>150 & POS<6601607' -o $i"_filtered_indels.vcf" --threads 8 $PF/raw_indels/$i"_raw_indels.vcf"; done

cd ..

printf "${COLOR}All done.${NC}\n"
