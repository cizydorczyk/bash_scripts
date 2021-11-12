# Set input and parameters:
assemblies_dir="/home/conrad/hinfluenzae/de_novo_assemblies/i16_nextpolish/i2_unpolished_assembly_links/"
fastq_dir="/home/conrad/hinfluenzae/fastq_files/i1_q5_trimmed_fastq/"
isolate_list="/home/conrad/hinfluenzae/hi_isolate_list.txt"

round=2
threads=8

count=0
for isolate in $(cat $isolate_list);
do
  assembly="${assemblies_dir}${isolate}_q5_unpolished_assembly.fasta"
  r1="${fastq_dir}${isolate}_1.fastq.gz"
  r2="${fastq_dir}${isolate}_2.fastq.gz"

  for ((i=1; i<=${round};i++))
  do
    echo "NOW STARTING ROUND ${i} for isolate ${isolate}..."
    # Step 1
    echo "Indexing original assembly for isolate ${isolate}..."
    bwa index ${assembly}
    echo "Aligning reads to assembly for isolate ${isolate}..."
    bwa mem -t ${threads} ${assembly} ${r1} ${r2} | samtools view --threads ${threads} -F 0x4 -b | samtools fixmate -m --threads ${threads} - - | samtools sort --threads ${threads} - | samtools markdup --threads ${threads} -r - ${isolate}.sgs.sort.bam
    # index bam and genome files
    echo "Indexing bam for isolate ${isolate}..."
    samtools index -@ ${threads} ${isolate}.sgs.sort.bam
    samtools faidx ${assembly}
    # polish genome file (round 1)
    echo "Polishing genome (task 1) for isolate ${isolate}..."
    python ~/software/NextPolish/lib/nextpolish1.py -g ${assembly} -t 1 -p ${threads} -s ${isolate}.sgs.sort.bam -ploidy 1 > ${isolate}.genome.polishtemp.fasta
    assembly=${isolate}.genome.polishtemp.fasta
    ## Step 2
    echo "Indexing task 1 output for isolate ${isolate}..."
    bwa index ${assembly}
    echo "Aligning reads to task 1 output assembly for isolate ${isolate}..."
    bwa mem -t ${threads} ${assembly} ${r1} ${r2} | samtools view --threads ${threads} -F 0x4 -b | samtools fixmate -m --threads ${threads} - - | samtools sort --threads ${threads} - | samtools markdup --threads ${threads} -r - ${isolate}.sgs.sort.bam
    # index bam and genome files
    echo "Indexing bam"
    samtools index -@ ${threads} ${isolate}.sgs.sort.bam
    samtools faidx ${assembly}
    # polish genome file (round 1)
    echo "Polishing genome (task 2)"
    python ~/software/NextPolish/lib/nextpolish1.py -g ${assembly} -t 2 -p ${threads} -s ${isolate}.sgs.sort.bam -ploidy 1 > ${isolate}_q5_polished_assembly.fasta
    assembly=${isolate}_q5_polished_assembly.fasta

    echo "FINISHED ROUND ${i}"
  done

  echo "Cleaning up......."

  rm *.bam *.bai *.amb *.ann *.bwt *.fai *.pac *.sa *temp.fasta

  echo "FINISHED PIPELINE."
done
# polished genome is ${isolate}_q5_polished_assembly.fasta
