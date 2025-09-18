#!/bin/bash

# 提取all.fimo.out.filt的第一列
cut -f1 all.fimo.out.filt > all.fimo_col1.txt
awk 'NR==FNR {a[$1]; next} $1 in a' all.fimo_col1.txt total/all.gene_id.txt > total/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' all.fimo_col1.txt a/a.gene_id.txt > a/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' all.fimo_col1.txt b/b.gene_id.txt > b/gene.TF_ids.txt

# 在第一行前添加标题行
sed -i '1iname\tDegree' total/gene.TF_ids.txt
sed -i '1iname\tDegree' a/gene.TF_ids.txt
sed -i '1iname\tDegree' b/gene.TF_ids.txt

# 清理临时文件
rm all.fimo_col1.txt
