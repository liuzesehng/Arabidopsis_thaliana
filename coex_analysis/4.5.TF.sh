#!/bin/bash

# 提取Ath_TF_list.txt的第二列
cut -f2 Ath_TF_list.txt > Ath_TF_list_col2.txt


awk 'NR==FNR {a[$1]; next} $1 in a' Ath_TF_list_col2.txt total/gene_id.tsv > total/gene.TF_ids.txt
grep -v -f total/gene.TF_ids.txt total/gene_id.tsv > total/gene.nonTF_ids.txt
# 在第一行前添加标题行
sed -i '1iname\tDegree' total/gene.TF_ids.txt


# 清理临时文件
rm Ath_TF_list_col2.txt
