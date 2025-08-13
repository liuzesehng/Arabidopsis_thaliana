#!/bin/bash
ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/salmon_index
tx2gene=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/mapping/tx2gene.tsv
for i in $(ls Normal/*/SR*gz | grep "SRR3465232_1")
do
    m=$(dirname $i)
    j=${i/_1/_2}
    # n=${m/_1}
    salmon quant -i $ref -g $tx2gene -l A -p 64 -1 $i -2 $j -o salmon3/$m 
done
#wait