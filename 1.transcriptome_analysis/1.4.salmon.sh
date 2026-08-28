#!/bin/bash
ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/salmon_index
tx2gene=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/mapping/tx2gene.tsv
for i in $(ls Abnormal/*/SR*gz)
do
    m=$(dirname $i)
    salmon quant -i $ref -g $tx2gene -l A -p 64 -r $i -o salmon3/$m
done

for i in $(ls Normal/*/SR*gz | grep -v "SRR3465232")
do
    m=$(dirname $i)
    salmon quant -i $ref -g $tx2gene -l A -p 64 -r $i -o salmon3/$m 
done
#wait
