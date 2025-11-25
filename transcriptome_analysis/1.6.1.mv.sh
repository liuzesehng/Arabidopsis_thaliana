#!/bin/bash
# 从biosample_result.txt中提取10C温度的GEO编号（不包含表头）
awk '/GEO:/{match($0,/GEO: ([A-Z]+[0-9]+)/,a); geo=a[1]} /growth temperature="10C"/{if(geo){print geo; geo=""}}' biosample_result.txt > 10C.csv
awk -F',' 'NR == 1' Unnormal.csv > RCA/Abnormal/10C/10C.csv
# 从Unnormal.csv中选取第30列匹配10C.csv第一列的所有行
awk -F',' 'FNR==NR{a[$1]; next} $30 in a' 10C.csv Unnormal.csv >> RCA/Abnormal/10C/10C.csv
awk -F',' '{print $1}' 10C.csv | grep -vFf - Unnormal.csv > RCA/Abnormal/16C/16C.csv
rm -rf 10C.csv

for i in $(ls RCA/Abnormal/*/*csv)
do
    m=$(dirname $i)
    if [ -d $m ]; then
        echo "dir exists"
    else
        mkdir $m
    fi
    for j in $(awk -F',' 'NR>1 {print $1}' $i; awk -F',' 'NR>1 {print $11}' $i)
    do
        #mv RCA/Abnormal/$j.RCA.transcript_count.txt $m
        mv RCA/Abnormal/$j.RCA.transcript.txt $m
    done
done
rm -rf RCA/Abnormal/10C/10C.csv RCA/Abnormal/16C/16C.csv
