(snippy-env) conrad@conrad-Precision-Tower-3620:~/ncbi_ecoli_genomes/ncbi_genomes_errFree_sim_reads$ for i in $(cat ~/ncbi_ecoli_genomes/ncbi-genomes-2019-06-07/in_silico_st131_typing/ncbi_st131_c1_isolates.txt); do ~/software/art_bin_MountRainier/art_illumina -ss MSv3 -p -i ~/ncbi_ecoli_genomes/ncbi-genomes-2019-06-07/ncbi_st73_131_1193_genomes/$i"_genomic.fna" -l 250 -m 400 -s 110 -f 75 -o $i -rs 42; cat $i"_errFree.sam" | grep -v ^@ | awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > $i"_1.fastq"; cat $i"_errFree.sam" | grep -v ^@ | awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > $i"_2.fastq"; rm $i".sam" $i"1.aln" $i"1.fq" $i"2.aln" $i"2.fq" $i"_errFree.sam"; done
