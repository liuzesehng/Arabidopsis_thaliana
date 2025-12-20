#!/bin/bash
# total Rca data
# for file in $(ls salmon3/Abnormal/*/quant.genes.sf)
# do
#     filename=$(awk -F'/' '{print $3}' <<< $file)
#     path=$(ls RCA/Abnormal/*/$filename*txt)
#     p=$(awk -v var="$filename" -F',' '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)
#     o=$(cat biosample_result.txt | grep -A 9 $p | grep "growth temperature=" | sed 's/.*\/growth temperature="\([^"]*\)".*/\1/')
#     q=$(cat biosample_result.txt | grep -A 9 $p | grep "accession number=" | sed 's/.*\/accession number="\([^"]*\)".*/\1/')

#     if [ "$q" = "991" ] && [ "$o" = "10C" ]; then
#         cut -f 1,4 $file > coex2/RCA.$o.tran.txt
#         awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/RCA.$o.tran.txt > temp && mv temp coex2/RCA.$o.tran.txt
#         # 对RCA.tran.txt进行排序（保留表头）
#         (head -n 1 coex2/RCA.$o.tran.txt; tail -n +2 coex2/RCA.$o.tran.txt | sort -k1,1) > temp && mv temp coex2/RCA.$o.tran.txt
#     elif [ "$q" = "991" ] && [ "$o" = "16C" ]; then
#         cut -f 1,4 $file > coex2/RCA.$o.tran.txt
#         awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/RCA.$o.tran.txt > temp && mv temp coex2/RCA.$o.tran.txt
#         # 对RCA.tran.txt进行排序（保留表头）
#         (head -n 1 coex2/RCA.$o.tran.txt; tail -n +2 coex2/RCA.$o.tran.txt | sort -k1,1) > temp && mv temp coex2/RCA.$o.tran.txt
#     else
#         cut -f 1,4 $file > coex2/$filename.$o.tran.txt
#         awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.$o.tran.txt > temp && mv temp coex2/$filename.$o.tran.txt
#         # 对两个文件进行排序（保留表头）然后join
#         (head -n 1 coex2/RCA.$o.tran.txt; tail -n +2 coex2/RCA.$o.tran.txt | sort -k1,1) > coex2/RCA.$o.tran.sorted.txt
#         (head -n 1 coex2/$filename.$o.tran.txt; tail -n +2 coex2/$filename.$o.tran.txt | sort -k1,1) > coex2/$filename.$o.tran.sorted.txt
#         join -t $'\t' -1 1 -2 1 coex2/RCA.$o.tran.sorted.txt coex2/$filename.$o.tran.sorted.txt > temp && mv temp coex2/RCA.$o.tran.txt
#         rm -rf coex2/$filename.$o.tran.txt coex2/RCA.$o.tran.sorted.txt coex2/$filename.$o.tran.sorted.txt
#     fi
# done

for file in $(ls salmon3/Normal/*/quant.genes.sf)
do
    filename=$(awk -F'/' '{print $3}' <<< $file)
    o=22C
    p=$(awk -v var="$filename" -F',' '($1 == var || $11 == var) { print $30 }' Normal.csv | uniq)
    q=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk -F'[()]' '{print $2}')

    if [ "$q" = "10008" ] && [ "$o" = "22C" ]; then
        cut -f 1,4 $file > coex2/RCA.$o.tran.txt
        awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/RCA.$o.tran.txt > temp && mv temp coex2/RCA.$o.tran.txt
        # 对RCA.tran.txt进行排序（保留表头）
        (head -n 1 coex2/RCA.$o.tran.txt; tail -n +2 coex2/RCA.$o.tran.txt | sort -k1,1) > temp && mv temp coex2/RCA.$o.tran.txt
    else
        cut -f 1,4 $file > coex2/$filename.$o.tran.txt
        awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.$o.tran.txt > temp && mv temp coex2/$filename.$o.tran.txt
        # 对两个文件进行排序（保留表头）然后join
        (head -n 1 coex2/RCA.$o.tran.txt; tail -n +2 coex2/RCA.$o.tran.txt | sort -k1,1) > coex2/RCA.$o.tran.sorted.txt
        (head -n 1 coex2/$filename.$o.tran.txt; tail -n +2 coex2/$filename.$o.tran.txt | sort -k1,1) > coex2/$filename.$o.tran.sorted.txt
        join -t $'\t' -1 1 -2 1 coex2/RCA.$o.tran.sorted.txt coex2/$filename.$o.tran.sorted.txt > temp && mv temp coex2/RCA.$o.tran.txt
        rm -rf coex2/$filename.$o.tran.txt coex2/RCA.$o.tran.sorted.txt coex2/$filename.$o.tran.sorted.txt
    fi
done

awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex2/RCA.10C.tran.txt | wc -l # 查看文件的列数
awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex2/RCA.16C.tran.txt | wc -l # 查看文件的列数
awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex2/RCA.22C.tran.txt | wc -l # 查看文件的列数

