#!/bin/bash
cat 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff | sed -n '12,12p' | awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > filtered_snp.txt
for i in $(awk '$1 <= 16573692' filtered_snp_frequencies.csv | cut -f 1)
do
    cat 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff | awk -v var="$i" '$1 == 2 && $2 == var' | awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' >> filtered_snp.txt
done