#!/bin/bash
ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/salmon_index
for i in $(ls Normal/SR*gz | grep "SRR3465232_1")
do
    m=${i/.fastq.gz}
    j=${i/_1/_2}
    n=${m/_1}
    salmon quant -i $ref -l A -p 64 -1 $i -2 $j -o salmon/$n 
done
#wait