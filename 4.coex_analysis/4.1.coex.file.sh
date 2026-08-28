#!/bin/bash
# total Rca data
for file in $(ls RCA/Abnormal/*/*.txt)
do
    filename=$(awk -F'/' '{print $3}' <<< $file)
    path=$(ls RCA3/Abnormal/*/*/*/$filename*txt)
    p=$(awk -v var="$filename" -F',' '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)
    o=$(cat biosample_result.txt | grep -A 9 $p | grep "growth temperature=" | sed 's/.*\/growth temperature="\([^"]*\)".*/\1/')
    q=$(cat biosample_result.txt | grep -A 9 $p | grep "accession number=" | sed 's/.*\/accession number="\([^"]*\)".*/\1/')

    if [ "$q" = "991" ]; then
        cut -f 1,4 $file > coex2/RCA.tran.txt
        awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/RCA.tran.txt > temp && mv temp coex2/RCA.tran.txt
        # 对RCA.tran.txt进行排序（保留表头）
        (head -n 1 coex2/RCA.tran.txt; tail -n +2 coex2/RCA.tran.txt | sort -k1,1) > temp && mv temp coex2/RCA.tran.txt
    else
        cut -f 1,4 $file > coex2/$filename.tran.txt
        awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.tran.txt > temp && mv temp coex2/$filename.tran.txt
        # 对两个文件进行排序（保留表头）然后join
        (head -n 1 coex2/RCA.tran.txt; tail -n +2 coex2/RCA.tran.txt | sort -k1,1) > coex2/RCA.tran.sorted.txt
        (head -n 1 coex2/$filename.tran.txt; tail -n +2 coex2/$filename.tran.txt | sort -k1,1) > coex2/$filename.tran.sorted.txt
        join -t $'\t' -1 1 -2 1 coex2/RCA.tran.sorted.txt coex2/$filename.tran.sorted.txt > temp && mv temp coex2/RCA.tran.txt
        rm -rf coex2/$filename.tran.txt coex2/RCA.tran.sorted.txt coex2/$filename.tran.sorted.txt
    fi
done

for file in $(ls salmon3/Normal/*/quant.genes.sf)
do
    filename=$(awk -F'/' '{print $3}' <<< $file)
    p=$(awk -v var="$filename" -F',' '($1 == var || $11 == var) { print $30 }' Normal.csv | uniq)
    q=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk -F'[()]' '{print $2}')
    cut -f 1,4 $file > coex2/$filename.tran.txt
    awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.tran.txt > temp && mv temp coex2/$filename.tran.txt
    # 对两个文件进行排序（保留表头）然后join
    (head -n 1 coex2/RCA.tran.txt; tail -n +2 coex2/RCA.tran.txt | sort -k1,1) > coex2/RCA.tran.sorted.txt
    (head -n 1 coex2/$filename.tran.txt; tail -n +2 coex2/$filename.tran.txt | sort -k1,1) > coex2/$filename.tran.sorted.txt
    join -t $'\t' -1 1 -2 1 coex2/RCA.tran.sorted.txt coex2/$filename.tran.sorted.txt > temp && mv temp coex2/RCA.tran.txt
    rm -rf coex2/$filename.tran.txt coex2/RCA.tran.sorted.txt coex2/$filename.tran.sorted.txt
done

awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex2/RCA.tran.txt | wc -l # 查看文件的列数


