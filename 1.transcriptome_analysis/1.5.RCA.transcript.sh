#!/bin/bash
mkdir -p RCA/Abnormal
mkdir -p RCA/Normal
for i in $(ls salmon3/Abnormal/*/quant.sf)
do
    n=$(basename $(dirname $i))
    awk 'NR==1 || $1=="NM_129531.3" || $1=="NM_179989.3" || $1=="NM_179990.1"' $i > RCA/Abnormal/$n.RCA.transcript.txt
    awk '$1=="AT2G39730"' salmon3/Abnormal/$n/quant.genes.sf >> RCA/Abnormal/$n.RCA.transcript.txt
done

for i in $(ls salmon3/Normal/*/quant.sf)
do
    n=$(basename $(dirname $i))
    awk 'NR==1 || $1=="NM_129531.3" || $1=="NM_179989.3" || $1=="NM_179990.1"' $i > RCA/Normal/$n.RCA.transcript.txt
    awk '$1=="AT2G39730"' salmon3/Normal/$n/quant.genes.sf >> RCA/Normal/$n.RCA.transcript.txt
done
