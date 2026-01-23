#!/bin/bash

# all.fimo.out.filt

awk 'NR==FNR {a[$1]; next} $1 in a' total/CytoscapeInput-nodes-brown.txt all.fimo.out.filt > total/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 10C/10C.CytoscapeInput-nodes-brown.txt all.fimo.out.filt > 10C/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 16C/16C.CytoscapeInput-nodes-red.txt all.fimo.out.filt > 16C/gene.TF_ids.txt
awk 'NR==FNR {a[$1]; next} $1 in a' 22C/22C.CytoscapeInput-nodes-green.txt all.fimo.out.filt > 22C/gene.TF_ids.txt

# 在第一行前添加标题行
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' total/gene.TF_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 10C/gene.TF_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 16C/gene.TF_ids.txt
sed -i '1imotif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence' 22C/gene.TF_ids.txt

# 清理临时文件
# rm all.fimo_col1.txt
