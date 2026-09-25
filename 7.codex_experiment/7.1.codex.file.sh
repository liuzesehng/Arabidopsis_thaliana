#!/bin/bash
# total Rca data
for file in $(ls rnadata/*/salmon/quant.genes.sf)
do
    filename=$(awk -F'/' '{print $2}' <<< $file)
    filename=${filename/T}
    if [ "$filename" = "12_2" ]; then
        cut -f 1,4 $file > coex_experiment/RCA.tran.txt
        awk -v new_col="$filename" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex_experiment/RCA.tran.txt > temp && mv temp coex_experiment/RCA.tran.txt
    else
        cut -f 1,4 $file > coex_experiment/$filename.tran.txt
        awk -v new_col="$filename" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex_experiment/$filename.tran.txt > temp && mv temp coex_experiment/$filename.tran.txt
        join -t $'\t' -1 1 -2 1 coex_experiment/RCA.tran.txt coex_experiment/$filename.tran.txt > temp && mv temp coex_experiment/RCA.tran.txt
        rm -rf coex_experiment/$filename.tran.txt
    fi
done

awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex_experiment/RCA.tran.txt | wc -l # 查看文件的列数


