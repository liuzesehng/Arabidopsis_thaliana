#!/bin/bash
ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/salmon_index
for i in $(ls Normal/SR*gz | grep -v "SRR3465232")
do
    m=${i/.fastq.gz}
    salmon quant -i $ref -l A -p 64 -r $i -o salmon/$m 
done
#wait