#!/bin/bash

while getopts 'p:f:t:i:' c
do
  case $c in
    p) PROJECTDIR=$OPTARG ;;
    f) INPUTFASTA=$OPTARG ;;
    t) TARGETSEQ=$OPTARG ;;
    i) ISOLATELIST=$OPTARG ;;
  esac
done

#cd ${PROJECTDIR}

# Get concatenated, rc masked reference:
#(base) conrad@conrad-Precision-Tower-3620:~/hinfluenzae/i43_hypermutator_testing/st393_public_ref/snp_calling/with_ref_iqtree_cfml_rcmasked_aln$ python ~/python_scripts/phd/F45_extract_seq_from_fasta_1.py --input_fasta rc_masked.clean.393.full.aln --output_fasta ../../tstv/reference.rc_masked.concat.fasta --target_seq Reference

echo "Obtaining concatenated reference sequence..."
CONCAT_REF=${PROJECTDIR}/concat_ref.fasta
python /home/conrad/python_scripts/phd/F45_extract_seq_from_fasta_1.py --input_fasta ${INPUTFASTA} --output_fasta ${CONCAT_REF} --target_seq ${TARGETSEQ}

# Split fasta files for VCF generation using snp-sites:
#(base) conrad@conrad-Precision-Tower-3620:~/hinfluenzae/i43_hypermutator_testing/st393_public_ref/tstv$ python ~/python_scripts/phd/F45_split_fasta_1.py --input_fasta ../snp_calling/with_ref_iqtree_cfml_rcmasked_aln/rc_masked.clean.393.full.aln --ref reference.rc_masked.concat.fasta --output_folder split_fasta_files/

echo "Splitting fasta files for VCF generation using snp-sites (must be installed)..."
SPLIT_FASTA_FILES_DIR=${PROJECTDIR}split_fasta_files/

mkdir -p ${SPLIT_FASTA_FILES_DIR}

python ~/python_scripts/phd/F45_split_fasta_1.py --input_fasta ${INPUTFASTA} --ref ${CONCAT_REF} --output_folder ${SPLIT_FASTA_FILES_DIR}

# Generate vcf files with only rc-masked SNPs using snp-sites:
#(base) conrad@conrad-Precision-Tower-3620:~/hinfluenzae/i43_hypermutator_testing/st393_public_ref/tstv$ for i in $(cat ~/hinfluenzae/st393_isolate_list.txt); do snp-sites -v -o isolate_vcf_files/$i".vcf" split_fasta_files/$i".fasta"; done

echo "Generating VCF files with only RC-masked SNPs using snp-sites..."
ISOLATE_VCF_FILES_DIR=${PROJECTDIR}isolate_vcf_files/

mkdir -p ${ISOLATE_VCF_FILES_DIR}

for i in $(cat ${ISOLATELIST}); do snp-sites -v -o ${ISOLATE_VCF_FILES_DIR}${i}.vcf ${SPLIT_FASTA_FILES_DIR}${i}.fasta; done

# Run VCFTools TsTv module:
#(base) conrad@conrad-Precision-Tower-3620:~/hinfluenzae/i43_hypermutator_testing/st393_public_ref/tstv$ for i in $(cat ~/hinfluenzae/st393_isolate_list.txt); do vcftools --vcf isolate_vcf_files/$i".vcf" --out vcftools_tstv_calc/$i --TsTv-summary; done

echo "Running VCFTools TsTv module..."
VCFTOOLS_TSTV_CALC_DIR=${PROJECTDIR}vcftools_tstv_calc/

mkdir -p ${VCFTOOLS_TSTV_CALC_DIR}

for i in $(cat ${ISOLATELIST}); do vcftools --vcf ${ISOLATE_VCF_FILES_DIR}${i}.vcf --out ${VCFTOOLS_TSTV_CALC_DIR}${i} --TsTv-summary; done

# Summarize VCFTools output:
#(base) conrad@conrad-Precision-Tower-3620:~/hinfluenzae/i43_hypermutator_testing/st393_public_ref/tstv$ for i in $(cat ~/hinfluenzae/st393_isolate_list.txt); do python ~/python_scripts/phd/i23_summarize_vcftools_tstv_1.py --input_log vcftools_tstv_calc/$i".log" --output_file tstv_summary.txt --isolate $i; done

echo "Summarizing VCFTools output..."
VCFTOOLS_OUTPUT_SUMMARY=${PROJECTDIR}tstv_summary.txt

for i in $(cat ${ISOLATELIST}); do python ~/python_scripts/phd/i23_summarize_vcftools_tstv_1.py --input_log ${VCFTOOLS_TSTV_CALC_DIR}${i}.log --output_file ${VCFTOOLS_OUTPUT_SUMMARY} --isolate ${i}; done

echo "Done."

