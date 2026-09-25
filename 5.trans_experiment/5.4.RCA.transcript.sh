#!/bin/bash
mkdir -p RCA/rnadata
for i in $(ls rnadata/*/salmon/*quant.sf)
do
    n=$(echo $i | cut -d'/' -f2)
    awk 'NR==1 || $1=="NM_129531.3" || $1=="NM_179989.3" || $1=="NM_179990.1"' $i > RCA/rnadata/$n.RCA.transcript.txt
    awk '$1=="AT2G39730"' rnadata/$n/salmon/quant.genes.sf >> RCA/rnadata/$n.RCA.transcript.txt
done
