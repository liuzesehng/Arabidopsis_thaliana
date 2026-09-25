#!/bin/bash
ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/salmon_index
tx2gene=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/mapping/tx2gene.tsv
for i in $(ls rnadata/*/T*_1.fit.fastq.gz)
do
    m=$(dirname $i)
    j=${i/_1.fit/_2.fit}
    salmon quant -i $ref -g $tx2gene -l A -p 64 -1 $i -2 $j -o $m/salmon
done
