#!/bin/bash
# 提取模块中的转录因子（TF）
awk 'NR==FNR {a[$1]; next} $2 in a' yellow_module_genes.txt ../coex2/Ath_TF_list.txt | cut -f2,3 | sort -u > TF_in_module.txt
sort -k2,2 TF_in_module.txt -o TF_in_module.txt


awk 'NR==FNR {a[$1]; next} $2 in a' negative_turquoise_module_genes.txt ../coex2/Ath_TF_list.txt | cut -f2,3 | sort -u > TF_negative_in_module.txt
sort -k2,2 TF_negative_in_module.txt -o TF_negative_in_module.txt


awk 'NR==FNR {a[$1]; next} $2 in a' other_module_genes_excluding_yellow_turquoise.txt ../coex2/Ath_TF_list.txt | cut -f2,3 | sort -u > TF_other_in_module.txt
sort -k2,2 TF_other_in_module.txt -o TF_other_in_module.txt

# 在第一行前添加标题行
sed -i '1iTF_ID\tTF_Family' TF_in_module.txt

sed -i '1iTF_ID\tTF_Family' TF_negative_in_module.txt

sed -i '1iTF_ID\tTF_Family' TF_other_in_module.txt

# all.fimo.out.filt 文件包含所有基因的TF结合位点信息
awk 'NR==FNR {a[$1]; next} $1 in a' yellow_module_genes.txt ../coex2/all.fimo.out.filt > gene.TF_ids.txt

awk 'NR==FNR {a[$1]; next} $1 in a' negative_turquoise_module_genes.txt ../coex2/all.fimo.out.filt > gene.TF_negative_ids.txt

awk 'NR==FNR {a[$1]; next} $1 in a' other_module_genes_excluding_yellow_turquoise.txt ../coex2/all.fimo.out.filt > gene.TF_other_ids.txt
# 在第一行前添加标题行
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' gene.TF_ids.txt

sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' gene.TF_negative_ids.txt

sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' gene.TF_other_ids.txt
# 清理临时文件
# rm all.fimo_col1.txt
