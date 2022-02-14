#!/bin/bash

##########################
### k5 bash script to get an annotated reference that includes contigs created using unmapped reads
### Written by: Conrad Izydorczyk
### Last updated: November 4, 2021
##########################
### Run in a conda environment with Unicycler and Bioawk installed.
### Rast annotation tools from PATRIC must also be available.
##########################

### Project directory must be directory with snp calling results (raw snippy output, snippy core, etc.)
### Isolate list is a text file with one isolate per line.
### ST is the sequence type being analyzed. Another idendifier may be used.
### MINCONTIGLEN = shortest contig length to keep among assembled, unmapped reads.
### REF = reference genome to which to append assembled contigs.
### WORKFLOW = rasttk workflow to use during annotation.
### GENUS = genus for which to keep contigs with blast hits to; genus must occur at least once in the top 10 hits

echo "Starting script."

while getopts 'p:i:s:m:r:n:w:g:' c
do
  case $c in
    p) PROJECTDIR=$OPTARG ;;
    i) ISOLATELIST=$OPTARG ;;
    s) ST=$OPTARG ;;
    m) MINCONTIGLEN=$OPTARG ;;
    r) REF=$OPTARG ;;
    n) SCINAME=$OPTARG ;;
    w) WORKFLOW=$OPTARG ;;
    g) GENUS=$OPTARG ;;
  esac
done

COMBFQDIR=${PROJECTDIR}combined-unmapped-reads/
UNIDIR=${COMBFQDIR}unicycler/
#
# mkdir -p ${COMBFQDIR}
#
# echo "Combining unmapped reads from all samples..."
#
COMBINEDR1=${COMBFQDIR}${ST}-comb-unmapped_1.fastq
COMBINEDR2=${COMBFQDIR}${ST}-comb-unmapped_2.fastq
#
# for i in $(cat ${ISOLATELIST})
# do
#     echo ${i}
#     UNMAPPEDR1=${PROJECTDIR}raw_snippy_output/${i}/snps.unmapped_R1.fq.gz
#     UNMAPPEDR2=${PROJECTDIR}raw_snippy_output/${i}/snps.unmapped_R2.fq.gz
#
#     zcat ${UNMAPPEDR1} >> ${COMBINEDR1}
#     zcat ${UNMAPPEDR2} >> ${COMBINEDR2}
# done
#
# gzip ${COMBINEDR1}
# gzip ${COMBINEDR2}
#
# echo "Obtaining read counts of combined forward & reverse read files..."
#
# fix_base_count() {
#     local counts=($(cat))
#     echo " ${counts[0]} $((${counts[1]} - ${counts[0]}))"
# }
#
# gzip -dc ${COMBINEDR1} \
#     | awk 'NR % 4 == 2' \
#     | wc -cl \
#     | fix_base_count
#
# gzip -dc ${COMBINEDR2} \
#     | awk 'NR % 4 == 2' \
#     | wc -cl \
#     | fix_base_count
#
# ### Assemble unmapped, combined reads using Unicycler:
# echo "Assembling unmapped reads using Unicycler..."
#
# unicycler -1 ${COMBINEDR1}.gz -2 ${COMBINEDR2}.gz -o ${UNIDIR} --depth_filter 0.25 --min_fasta_len 100 --min_polish_size 100 -t 8
#
UNIASSEMBLY=${UNIDIR}assembly.fasta
FILTASSEMBLY=${UNIDIR}${ST}-filtered-assembly.fasta
#
# # ### Filter assembled contigs based on length (move to python script? filtering pre-blast might be good though to remove garbage...):
# # echo "Filtering contigs shorter than minimum length provided..."
#
# bioawk -v len=${MINCONTIGLEN} -c fastx '{ if(length($seq) > len) { print ">"$name; print $seq }}' ${UNIASSEMBLY} > ${FILTASSEMBLY}

### Blast assembled contigs & remove non-target genus contigs
RAWBLASTOUTPUT=${UNIDIR}${ST}-raw-blast-output.xml

echo "Running blastn on contigs..."
# blastn -db nt -query ${FILTASSEMBLY} -outfmt "6 qseqid sseqid sacc bitscore evalue pident length staxids sscinames" -out ${RAWBLASTOUTPUT} -evalue 0.00001 -num_threads 8
blastn -db /media/conrad/Secondary_HDD/blastdb/nt -query ${FILTASSEMBLY} -outfmt 5 -evalue 0.00001 -out ${RAWBLASTOUTPUT} -num_threads 8 -max_target_seqs 10


### Parse blast output, keeping only contigs that have GENUS in at least one of their top 10 hits:
CONCATREF=${PROJECTDIR}${ST}-ref-with-unmapped-contigs.fasta

python /home/conrad/python_scripts/phd/m5-parse-blast-xml.py --blast_file ${RAWBLASTOUTPUT} --genus ${GENUS} --blast_contigs ${FILTASSEMBLY} --reference ${REF} --output_contigs ${CONCATREF}

### Concatenate ST-specific reference and filtered contigs from unmapped reads:
# CONCATREF=${PROJECTDIR}${ST}-ref-with-unmapped-contigs.fasta
# TEMPCONCATREF=${PROJECTDIR}${ST}-temp.fasta

# cat ${REF} ${FILTASSEMBLY} > ${TEMPCONCATREF}

# # renumber contigs in concatref b/c the added ones will repeat 1, 2, 3, ...
# awk '/^>/{print ">" ++i; next}{print}' < ${TEMPCONCATREF} > ${CONCATREF}
# rm ${TEMPCONCATREF}

### Annotate concatenated reference with rasttk:
echo "Annotating concatenated reference genome with rasttk..."

RASTDIR=${PROJECTDIR}rast-annotation/
RAWGTO=${RASTDIR}${ST}-ref-with-unmapped-contigs.gto
ANNGTO=${RASTDIR}${ST}-ref-with-unmapped-contigs-annotated.gto
ANNGFF=${RASTDIR}${ST}.gff
ANNGBK=${RASTDIR}${ST}.gbk

mkdir -p ${RASTDIR}

echo "Creating gto object..."
rast-create-genome --scientific-name ${SCINAME} --genetic-code 11 --domain Bacteria --contigs ${CONCATREF} > ${RAWGTO}

echo "Annotating gto object..."
rast-process-genome --workflow ${WORKFLOW} < ${RAWGTO} > ${ANNGTO}

echo "Exporting gff and gbk files..."
rast-export-genome gff < ${ANNGTO} > ${ANNGFF}
rast-export-genome genbank < ${ANNGTO} > ${ANNGBK}

echo "Done annotating."
