#!/bin/bash

# Programs must be in path:
# samtools
# bwa mem

# Program paths to set:
PICARD=/home/conrad/software/picard.jar
PLATANUS_B=/home/conrad/software/Platanus_B_v1.3.2/platanus_b

# File paths to set:
prodigal_training_file=/home/conrad/hinfluenzae/virulence_factors/F56_iga/alignments/rdkw20-prodigal.trn

# Set input flags
# -r (reference) must already be indexed by bwa
while getopts d:i:t:1:2:r: flag
do
    case "${flag}" in
        d) wd=${OPTARG};;
        i) isolate=${OPTARG};;
        t) threads=${OPTARG};;
	1) fastq1=${OPTARG};;
	2) fastq2=${OPTARG};;
	r) reference=${OPTARG};;
    esac
done

### make working directories if do not exist:
mkdir -p $wd

fastq_files_dir=${wd}fastq_files
assemblies_dir=${wd}assemblies

mkdir -p ${fastq_files_dir}
mkdir -p ${assemblies_dir}

### Get fastq files from alignments ###

## cd to fastq_files_dir to run alignments/generate fastq files
cd ${fastq_files_dir}

## run bwa
raw_aln_sam=${isolate}_raw_aln.sam

echo "Running bwa mem..."
bwa mem -t ${threads} ${reference} ${fastq1} ${fastq2} > ${raw_aln_sam}
echo "Done."

## extract mapped reads using samtools
sorted_bam=${isolate}_sorted.bam
unmapped_but_mate_mapped_bam=${isolate}_unmapped_but_mate_mapped.bam
mapped_but_mate_unmapped_bam=${isolate}_mapped_but_mate_unmapped.bam
both_reads_mapped_bam=${isolate}_both_reads_mapped.bam
all_mapped_sorted_reads_bam=${isolate}_all_mapped_sorted_reads.bam

# sort & convert to bam
samtools view -@ ${threads} -O BAM ${raw_aln_sam} | samtools sort -@ ${threads} - > ${sorted_bam}

# get mapped reads
samtools view -@ ${threads} -b -f 4 -F 8 ${sorted_bam} > ${unmapped_but_mate_mapped_bam}
samtools view -@ ${threads} -b -f 8 -F 4 ${sorted_bam} > ${mapped_but_mate_unmapped_bam}
samtools view -@ ${threads} -b -F 12 ${sorted_bam} > ${both_reads_mapped_bam}

samtools merge - ${unmapped_but_mate_mapped_bam} ${mapped_but_mate_unmapped_bam} ${both_reads_mapped_bam} | samtools sort -@ ${threads} - > ${all_mapped_sorted_reads_bam}

# extract reads from bam
mapped_f1=${fastq_files_dir}/${isolate}_mapped_1.fastq
mapped_f2=${fastq_files_dir}/${isolate}_mapped_2.fastq

java -Xmx20g -jar ${PICARD} SamToFastq I=${all_mapped_sorted_reads_bam} F=${mapped_f1} F2=${mapped_f2}

## clean up alignment files
rm ${sorted_bam} ${unmapped_but_mate_mapped_bam} ${mapped_but_mate_unmapped_bam} ${both_reads_mapped_bam} ${all_mapped_sorted_reads_bam} ${raw_aln_sam}

### Assemble reads ###

## cd to assemblies directory
cd ${assemblies_dir}

## mkdir for skesa and platanus_b output; unicycler creates its own output directory
platanus_b_output_dir=${assemblies_dir}/platanus_b_output
unicycler_output_dir=${assemblies_dir}/${isolate}_unicycler_output

mkdir -p ${platanus_b_output_dir}

## run unicycler
echo "Running Unicycler..."
unicycler -1 ${mapped_f1} -2 ${mapped_f2} --min_fasta_length 200 -t ${threads} -o ${unicycler_output_dir}
echo "Done running Unicycler."

# rename unicycler assembly
mv ${unicycler_output_dir}/assembly.fasta ${assemblies_dir}/${isolate}_unicycler_assembly.fasta

## run skesa
echo "Running skesa..."
skesa --fastq ${mapped_f1},${mapped_f2} --cores ${threads} --memory 32 > ${assemblies_dir}/${isolate}_skesa_assembly.fasta

## run platanus_b
platanus_b_raw_assembly=${platanus_b_output_dir}/${isolate}_platanus_b_raw_assembly
platanus_b_final_assembly=${platanus_b_output_dir}/${isolate}_platanus_b_assembly

echo "Running platanus_b..."
${PLATANUS_B} assemble -f ${mapped_f1} ${mapped_f2} -t ${threads} -m 32 -o ${platanus_b_raw_assembly}
${PLATANUS_B} iterate -c ${platanus_b_raw_assembly}_contig.fa -IP1 ${mapped_f1} ${mapped_f2} -i 10 -t ${threads} -tmp ${platanus_b_output_dir}

# mv platanus_b assembly
mv ${assemblies_dir}/out_iterativeAssembly.fa ${isolate}_platanus_b_assembly.fasta

## clean up assembly files
rm -r ${unicycler_output_dir} ${assemblies_dir}/out_iterateIntermediateResults ${platanus_b_output_dir}

echo "Done assemblies."

### Predict ORFs ###

echo "Running ORF prediction using Prodigal..."

## mkdir for orf predictions
orfs_dir=${wd}/prodigal_orfs

mkdir -p ${orfs_dir}

## run prodigal on assemblies
# NOTE: training file must be created manually & specified here

prodigal_output_prefix=${orfs_dir}/${isolate}

prodigal -i ${assemblies_dir}/${isolate}_unicycler_assembly.fasta -t ${prodigal_training_file} -o ${prodigal_output_prefix}_unicycler_orfs.gff -p single -f gff -d ${prodigal_output_prefix}_unicycler_orfs.fna -a ${prodigal_output_prefix}_unicycler_orfs.faa

prodigal -i ${assemblies_dir}/${isolate}_skesa_assembly.fasta -t ${prodigal_training_file} -o ${prodigal_output_prefix}_skesa_orfs.gff -p single -f gff -d ${prodigal_output_prefix}_skesa_orfs.fna -a ${prodigal_output_prefix}_skesa_orfs.faa

prodigal -i ${assemblies_dir}/${isolate}_platanus_b_assembly.fasta -t ${prodigal_training_file} -o ${prodigal_output_prefix}_platanus_b_orfs.gff -p single -f gff -d ${prodigal_output_prefix}_platanus_b_orfs.fna -a ${prodigal_output_prefix}_platanus_b_orfs.faa

echo "ORF prediction done."

