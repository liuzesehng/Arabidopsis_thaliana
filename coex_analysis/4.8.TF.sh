#!/bin/bash
# 提取模块中的转录因子（TF）
awk 'NR==FNR {a[$1]; next} $2 in a' total/yellow_module_genes.txt Ath_TF_list.txt > total/TF_in_module.txt
sort -k2,2 total/TF_in_module.txt -o total/TF_in_module.txt
awk 'NR==FNR {a[$1]; next} $2 in a' 10C/10C.blue_module_genes.txt Ath_TF_list.txt > 10C/TF_in_module.txt
sort -k2,2 10C/TF_in_module.txt -o 10C/TF_in_module.txt
awk 'NR==FNR {a[$1]; next} $2 in a' 16C/16C.blue_module_genes.txt Ath_TF_list.txt > 16C/TF_in_module.txt
sort -k2,2 16C/TF_in_module.txt -o 16C/TF_in_module.txt
awk 'NR==FNR {a[$1]; next} $2 in a' 22C/22C.yellow_module_genes.txt Ath_TF_list.txt > 22C/TF_in_module.txt
sort -k2,2 22C/TF_in_module.txt -o 22C/TF_in_module.txt

awk 'NR==FNR {a[$1]; next} $2 in a' total/negative_grey_module_genes.txt Ath_TF_list.txt > total/TF_negative_in_module.txt
sort -k2,2 total/TF_negative_in_module.txt -o total/TF_negative_in_module.txt
awk 'NR==FNR {a[$1]; next} $2 in a' 10C/10C.negative_grey_module_genes.txt Ath_TF_list.txt > 10C/TF_negative_in_module.txt
sort -k2,2 10C/TF_negative_in_module.txt -o 10C/TF_negative_in_module.txt
awk 'NR==FNR {a[$1]; next} $2 in a' 16C/16C.negative_pink_module_genes.txt Ath_TF_list.txt > 16C/TF_negative_in_module.txt
sort -k2,2 16C/TF_negative_in_module.txt -o 16C/TF_negative_in_module.txt
awk 'NR==FNR {a[$1]; next} $2 in a' 22C/22C.negative_black_module_genes.txt Ath_TF_list.txt > 22C/TF_negative_in_module.txt
sort -k2,2 22C/TF_negative_in_module.txt -o 22C/TF_negative_in_module.txt

# 在第一行前添加标题行
sed -i '1iTF_ID\tGene_ID\tTF_Family' total/TF_in_module.txt
sed -i '1iTF_ID\tGene_ID\tTF_Family' 10C/TF_in_module.txt
sed -i '1iTF_ID\tGene_ID\tTF_Family' 16C/TF_in_module.txt
sed -i '1iTF_ID\tGene_ID\tTF_Family' 22C/TF_in_module.txt

sed -i '1iGene_ID\tTF_ID\tTF_Family' total/TF_negative_in_module.txt
sed -i '1iGene_ID\tTF_ID\tTF_Family' 10C/TF_negative_in_module.txt
sed -i '1iGene_ID\tTF_ID\tTF_Family' 16C/TF_negative_in_module.txt
sed -i '1iGene_ID\tTF_ID\tTF_Family' 22C/TF_negative_in_module.txt

# all.fimo.out.filt 文件包含所有基因的TF结合位点信息
awk 'NR==FNR {a[$1]; next} $1 in a' total/yellow_module_genes.txt all.fimo.out.filt > total/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 10C/10C.blue_module_genes.txt all.fimo.out.filt > 10C/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 16C/16C.blue_module_genes.txt all.fimo.out.filt > 16C/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 22C/22C.yellow_module_genes.txt all.fimo.out.filt > 22C/gene.TF_ids.txt

awk 'NR==FNR {a[$1]; next} $1 in a' total/negative_grey_module_genes.txt all.fimo.out.filt > total/gene.TF_negative_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 10C/10C.negative_grey_module_genes.txt all.fimo.out.filt > 10C/gene.TF_negative_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 16C/16C.negative_pink_module_genes.txt all.fimo.out.filt > 16C/gene.TF_negative_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 22C/22C.negative_black_module_genes.txt all.fimo.out.filt > 22C/gene.TF_negative_ids.txt
# 在第一行前添加标题行
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' total/gene.TF_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 10C/gene.TF_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 16C/gene.TF_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 22C/gene.TF_ids.txt

sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' total/gene.TF_negative_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 10C/gene.TF_negative_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 16C/gene.TF_negative_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 22C/gene.TF_negative_ids.txt
# 清理临时文件
# rm all.fimo_col1.txt
