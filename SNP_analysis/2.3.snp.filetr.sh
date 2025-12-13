#!/bin/bash
cat 1001_SNP_MATRIX/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff | sed -n '12,12p' | awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > rca_snpeff.txt
for i in $(cut -f 2 filtered_snp_frequencies.csv)
do
    cat 1001_SNP_MATRIX/1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff | awk -v var="$i" '$1 == 2 && $2 == var' | awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' >> rca_snpeff.txt
done
