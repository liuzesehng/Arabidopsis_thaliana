#!/bin/bash
#conda activate meme
mkdir Ath_TF_binding_motifs
process_motifs() {
	local file=$1
	local dir=$(basename $file)
	local dir=${dir/.meme/_fimo_out}
	fimo --max-stored-scores 500000 --oc Ath_TF_binding_motifs/$dir $file rca.fasta
}
export -f process_motifs

ls Ath_TF_binding_motifs_individual/*.meme | parallel -j 64 process_motifs {}

# 合并所有fimo，过滤空行和p值小于0.00001的记录
cat Ath_TF_binding_motifs/*fimo_out/fimo.tsv | awk 'NF>0 && $8<0.00001' > all.fimo.out.filt
