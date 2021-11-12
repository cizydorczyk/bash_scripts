#!/bin/bash

# Uncomment lines in loop to run entire script
# Comment out lines to run only specific steps/in a stepwise manner
WORKDIR=/home/conrad/seq-lib-testing/de-novo-assemblies/m7-depth-length-plots/
FASTQDIR=/home/conrad/seq-lib-testing/fastq-files/m2-trimmed-reads/
ISOLATELIST=/home/conrad/seq-lib-testing/samples-list.txt
ASSEMBLYDIR=/home/conrad/seq-lib-testing/de-novo-assemblies/m5-unicycler-assemblies/unicycler-assemblies/

for i in $(cat ${ISOLATELIST})
do
	#ln -s ${ASSEMBLYDIR}${i}/${i}-assembly.fasta ${WORKDIR}${i}-assembly.fasta	
	#bwa index ${WORKDIR}${i}-assembly.fasta
	#bwa mem -t 8 ${WORKDIR}${i}-assembly.fasta ${FASTQDIR}${i}-1.fastq.gz ${FASTQDIR}${i}-2.fastq.gz > ${WORKDIR}${i}.sam
	#samtools view -@ 8 -O BAM -o ${WORKDIR}${i}.bam ${WORKDIR}${i}.sam
	#rm ${WORKDIR}${i}.sam
	#samtools sort -@ 8 ${WORKDIR}${i}.bam > ${WORKDIR}${i}.sorted.bam
	#rm ${WORKDIR}${i}.bam
	/home/conrad/software/bamtools/src/bamtools split -in ${WORKDIR}${i}.sorted.bam -reference
	for j in ${WORKDIR}${i}.sorted.REF_[0-9]*.bam; do echo ${j} >> ${WORKDIR}${i}.contig-depths.txt; samtools depth -a ${j} | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' >> ${WORKDIR}${i}.contig-depths.txt; done
	rm *REF_*
done
