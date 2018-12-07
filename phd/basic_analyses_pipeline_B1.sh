#!/bin/bash

## Conda must be installed & in path!
## Python must be in path! (Python 3)
## MLST Pipeline (v3) requires Torsten Seemann's 'mlst' and 'blastn' to be in path!
##	MLST Pipeline (v3) requires a directory with a fasta file with all allele sequences per gene
##	(single directory containing all fasta files)! Set path to this directory below!
## Fastq files (raw) must be unzipped!


## Set tool directories:
fastqc=/home/conrad/software/FastQC/fastqc
## Path to trimmomatic (installed with Anaconda in this case in default env):
trimmomatic=trimmomatic
## Must know which adapters to trim (run fastqc separately prior to running this pipeline):
trimmomatic_adapters_file=/home/conrad/software/Trimmomatic-0.38/adapters/NexteraPE-PE.fa
## Name of Shovill conda environment:
shovill_env=shovill-env
## Path to Shovill (installed with Anaconda in this case):
shovill=shovill
## Path to multiqc (installed with pip in this case):
multiqc=multiqc
## Path to kraken (v2 in this case):
kraken=/home/conrad/software/kraken2-2.0.7-beta/kraken2
## Path to kraken db to use:
kraken_db=/home/conrad/software/kraken2-2.0.7-beta/minikraken_8GB_v1
## Path to QUAST (installed with Anaconda in this case)
quast=quast
## Path to QUAST conda environment:
quast_env=quast-env
## Path to python script to estimate coverages:
coverage_py=/home/conrad/python_scripts/phd/get-total-bases-from-fastq-1.py
## Path to mlst pipeline script (currently v3)**SEE NOTE AT TOP OF SCRIPT:
mlst_pipeline_py=/home/conrad/python_scripts/phd/mlst-pipeline-3.py
## Path to directory with mlst allele fasta files:
mlst_alleles_fasta_dir=/home/conrad/hinfluenzae/basic_analysis_pipeline_test
## Comma separated list of mlst genes (same as names of fasta files with alleles):
mlst_alleles=adk,atpg,frdb,fuck,mdh,pgi,reca
## MLST scheme to use with Torsten Seemann's 'mlst' program:
mlst_scheme=hinfluenzae
## Set length of concatenated mlst genes (for QC of MLSA alignment). Shorter/longer isolate sequences will be removed from MLSA phylogenetic tree:
mlst_length=3057
## Set path to MLSA filtering script:
mlsa_filter_py=/home/conrad/python_scripts/phd/check-mlsa-lengths-1.py
## Set path to IQ-Tree:
iqtree=/home/conrad/software/iqtree-1.6.8-Linux/bin/iqtree

## Input variables:
## 1 = isolate list (full path); text file with one isolate per line
## 2 = project directory (full path); where all output will be generated
## 3 = fastq directory (full path); fastq files must be named like 'isolate-number_R1.fastq'
## 4 = estimated genome length of species (to estimate coverage)

## Make top level project directories:
cd $2
mkdir fastq_files fastqc de_novo_assemblies kraken mlst snp_calling phylogenetic_analyses

## Run trimmomatic:

mkdir $2/fastq_files/trimmed_fastq
mkdir $2/fastq_files/trimmed_fastq/paired
mkdir $2/fastq_files/trimmed_fastq/unpaired
mkdir $2/fastq_files/trimmed_fastq/log_files
mkdir $2/fastq_files/trimmed_fastq/summary_files

for i in $(cat $1)
do 
	$trimmomatic PE -threads 8 -trimlog $2/fastq_files/trimmed_fastq/log_files/$i"_log.txt" -summary $2/fastq_files/trimmed_fastq/summary_files/$i"_summary.txt" $3/$i"_R1.fastq" $3/$i"_R2.fastq" $2/fastq_files/trimmed_fastq/paired/$i"_paired_R1.fastq" $2/fastq_files/trimmed_fastq/unpaired/$i"_unpaired_R1.fastq" $2/fastq_files/trimmed_fastq/paired/$i"_paired_R2.fastq" $2/fastq_files/trimmed_fastq/unpaired/$i"_unpaired_R2.fastq" CROP:300 SLIDINGWINDOW:4:5 ILLUMINACLIP:$trimmomatic_adapters_file:2:30:10 MINLEN:50
done

## Run fastqc/multiqc on trimmed files:
mkdir $2/fastqc/trimmed_fastqc/

for i in $(cat $1)
do
	$fastqc $2/fastq_files/trimmed_fastq/paired/$i"_paired_R1.fastq"
	$fastqc $2/fastq_files/trimmed_fastq/paired/$i"_paired_R2.fastq"
	mv $2/fastq_files/trimmed_fastq/paired/*fastqc* $2/fastqc/trimmed_fastqc/
done

multiqc $2/fastqc/trimmed_fastqc/
mv $2/multiqc_report.html $2/multiqc_data/
mv $2/multiqc_data/ $2/fastqc/

## Run Kraken on trimmed fastq files:
mkdir $2/kraken/raw_kraken_output/
mkdir $2/kraken/kraken_reports/

for i in $(cat $1)
do
	$kraken --db $kraken_db --threads 8 --output $2/kraken/raw_kraken_output/$i"_kraken_output.txt" --report $2/kraken/kraken_reports/$i"_kraken_report.txt" --paired --use-names $2/fastq_files/trimmed_fastq/paired/$i"_paired_R1.fastq" $2/fastq_files/trimmed_fastq/paired/$i"_paired_R2.fastq"
done

## Assemble genomes with shovill:
mkdir $2/de_novo_assemblies/shovill_assemblies
source activate $shovill_env

for i in $(cat $1)
do
	$shovill --minlen 200 --cpus 8 --assembler spades --opts "--careful" --outdir $2/de_novo_assemblies/shovill_assemblies/$i/ --R1 $2/fastq_files/trimmed_fastq/paired/$i"_paired_R1.fastq" --R2 $2/fastq_files/trimmed_fastq/paired/$i"_paired_R2.fastq"
done
source deactivate

## Rename assembled contigs to include isolate numbers & run QUAST:
for i in $(cat $1)
do
	mv $2/de_novo_assemblies/shovill_assemblies/$i/contigs.fa $2/de_novo_assemblies/shovill_assemblies/$i/$i"_contigs.fa"
done

mkdir $2/de_novo_assemblies/quast
mkdir $2/de_novo_assemblies/quast/input_data
mkdir $2/de_novo_assemblies/quast/quast_output

for i in $(cat $1)
do
	cp $2/de_novo_assemblies/shovill_assemblies/$i/$i"_contigs.fa" $2/de_novo_assemblies/quast/input_data
done

source activate $quast_env

$quast -o $2/de_novo_assemblies/quast/quast_output/ -t 8 -m 200 $2/de_novo_assemblies/quast/input_data/*.fa
source deactivate

## Estimate isolate coverages:

for i in $(cat $1)
do
	python $coverage_py $2/fastq_files/trimmed_fastq/paired/$i"_paired_R1.fastq" $2/fastq_files/trimmed_fastq/paired/$i"_paired_R2.fastq" $4 $2/de_novo_assemblies/isolate_coverages.txt
done

## Run MLST+MLSA:
mv $2/de_novo_assemblies/quast/input_data $2/mlst/input_data
mkdir $2/mlst/mlst_pipeline_output

python $mlst_pipeline_py --mlsa --mlst_genes $mlst_alleles --mlst_alleles_directory $mlst_alleles_fasta_dir --scheme $mlst_scheme $2/mlst/input_data $2/mlst/mlst_pipeline_output

## Check MLSA alignment to make sure all isolates have proper alignments; deviate sequences are discarded (see python script for details):
mkdir $2/mlst/filtered_mlsa

python $mlsa_filter_py $2/mlst/mlst_pipeline_output/mlsa_alignment.fasta $2/mlst/mlst_pipeline_output/TS_mlst_results.txt $2/mlst/filtered_mlsa/filtered_mlsa.fasta $mlst_length $2/mlst/filtered_mlsa/filter_mlsa_log.txt $2/mlst/filtered_mlsa/removed_alignments.fasta

## Build ML tree using IQ-Tree; will use model select from IQ-Tree:
mkdir $2/phylogenetic_analyses/mlst_iqtree

$iqtree -nt AUTO -bb 10000 -m MFP -seed 142536 -pre $2/phylogenetic_analyses/mlst_iqtree/iq_tree_mlst -s $2/mlst/filtered_mlsa/filtered_mlsa.fasta



































