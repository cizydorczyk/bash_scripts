#!/bin/bash

# Set tool paths:

BWA=/home/conrad/Software/bwa-0.7.15/bwa
SAMTOOLS=/home/conrad/Software/samtools-1.3.1/samtools
BCFTOOLS=/home/conrad/Software/bcftools-1.3.1/bcftools

# Set project folder:

PF=/home/conrad/Data/primary_project_3/reference_alignments/H18_lesb58_indels

# Set reads folder:

READS=/home/conrad/Data/primary_project_3/fastq_files/G37_good_reads_rarefied

# Set printed-to-screen script text color:
COLOR='\033[1;36m'
NC='\033[0m'

# Input variables:

# 1 = reference genome (fasta file)
# 2 = isolate list (text file, one isolate number per line)

cd $PF

mkdir reference

cp $1 $PF/reference/reference.fasta

cd reference

echo "Reference genome is: "$1 > reference_READ_ME.txt

# Index reference genome using bwa index:

printf "${COLOR}Building bwa index of reference...${NC}\n"
$BWA index reference.fasta

# Index reference using samtools faidx here so as to not have to go back and do it after sorting/indexing bam files later:
printf "${COLOR}Building samtools faidx of reference...${NC}\n"
$SAMTOOLS faidx reference.fasta

cd ..

mkdir sorted_indexed_bam_alignments

cd sorted_indexed_bam_alignments

# Align reads to reference using bwa, convert output sam to bam, and sort the bam:
printf "${COLOR}Aligning reads to reference and sorting bam file...${NC}\n"
for i in $(cat $2); do $BWA mem -t 8 $PF/reference/reference.fasta $READS/$i"_"*"R1.fastq" $READS/$i"_"*"R2.fastq" | $SAMTOOLS view -bhuS - | $SAMTOOLS sort -@ 8 -O bam -o $i"_sorted_aln.bam" -; done

# Index bam files using samtools index:

printf "${COLOR}Indexing bam files...${NC}\n"
for i in *.bam; do $SAMTOOLS index $i; done

# cd ..

# mkdir raw_vcf_files

# cd raw_vcf_files

# # Call variants using samtools mpileup followed by bcftools call:

# printf "${COLOR}Calling variants...${NC}\n"
# for i in $(cat $2); do $SAMTOOLS mpileup -ugf $PF/reference/reference.fasta $PF/sorted_indexed_bam_alignments/$i"_sorted_aln.bam" | $BCFTOOLS call -m -v --ploidy 1 --threads 8 -o $i"_raw_variants.vcf" -; done
