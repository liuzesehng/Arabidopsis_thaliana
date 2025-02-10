#!/bin/bash
mkdir -p RCA/Normal
for i in $(ls salmon/Normal/*/quant.sf)
do
    n=$(basename $(dirname $i))
    awk 'NR==1 || $1=="NM_129531.3" || $1=="NM_179989.3" || $1=="NM_179990.1"' $i > RCA/Normal/$n.RCA.transcript.txt
done