#!/bin/bash
# coex_group
# 添加标志变量，用于控制per2值的显示
per2_shown=false

for file in $(ls salmon3/Abnormal/*/quant.genes.sf)
do
    filename=$(awk -F'/' '{print $3}' <<< $file)
    path=$(ls RCA3/Abnormal/*/*/*/$filename*txt)
    p=$(awk -v var="$filename" -F',' '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)
    o=$(cat biosample_result.txt | grep -A 9 $p | grep "growth temperature=" | sed 's/.*\/growth temperature="\([^"]*\)".*/\1/')
    q=$(cat biosample_result.txt | grep -A 9 $p | grep "accession number=" | sed 's/.*\/accession number="\([^"]*\)".*/\1/')
    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR==5 {print $4}' $path)

    # 计算第四列第二行的值
    value2=$(awk 'NR==2 {print $4}' $path)

    # 计算第四列第三行的值
    value3=$(awk 'NR==3 {print $4}' $path)

    # 计算第四列第四行的值
    value4=$(awk 'NR==4 {print $4}' $path)

    # 计算第四列第三行到第四行的值的总和
    sum2=$(awk 'NR>=3 && NR<=4 {sum += $4} END {print sum}' $path)

    # 调试输出
    echo "Processing file: $path"
    echo "Values: sum=$sum, value2=$value2, value3=$value3, value4=$value4"

    # 确保 value2, value3, value4 和 sum 不是空的
    if [ -n "$value2" ] && [ -n "$value3" ] && [ -n "$value4" ] && [ -n "$sum" ]; then
        # 确保所有变量转换为普通浮点数格式
        sum=$(echo $sum | awk '{printf "%.10f\n", $1}')
        sum2=$(echo $sum2 | awk '{printf "%.10f\n", $1}')
        value2=$(echo $value2 | awk '{printf "%.10f\n", $1}')
        value3=$(echo $value3 | awk '{printf "%.10f\n", $1}')
        value4=$(echo $value4 | awk '{printf "%.10f\n", $1}')

        # 使用 bc 计算
        persum2=$(echo "scale=4; ($sum2/$sum)*100" | bc)
        per2=$(echo "scale=4; ($value2/$sum)*100" | bc)
        per3=$(echo "scale=4; ($value3/$sum)*100" | bc)
        per4=$(echo "scale=4; ($value4/$sum)*100" | bc)

        # 添加前导零
        persum2=$(printf "%.4f" $persum2)
        per2=$(printf "%.4f" $per2)
        per3=$(printf "%.4f" $per3)
        per4=$(printf "%.4f" $per4)
    else
        echo "Warning: Missing values in $path"
        per2=0
        per3=0
        per4=0
    fi

    if [ $(echo "$per2 < 50" | bc) -eq 1 ]; then
        if [ "$q" = "991" ]; then
            cut -f 1,4 $file > coex2/RCA.b.tran.txt
            awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/RCA.b.tran.txt > temp && mv temp coex2/RCA.b.tran.txt
            # 对RCA.tran.txt进行排序（保留表头）
            (head -n 1 coex2/RCA.b.tran.txt; tail -n +2 coex2/RCA.b.tran.txt | sort -k1,1) > temp && mv temp coex2/RCA.b.tran.txt
        else
            cut -f 1,4 $file > coex2/$filename.b.tran.txt
            awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.b.tran.txt > temp && mv temp coex2/$filename.b.tran.txt
            # 对两个文件进行排序（保留表头）然后join
            (head -n 1 coex2/RCA.b.tran.txt; tail -n +2 coex2/RCA.b.tran.txt | sort -k1,1) > coex2/RCA.b.tran.sorted.txt
            (head -n 1 coex2/$filename.b.tran.txt; tail -n +2 coex2/$filename.b.tran.txt | sort -k1,1) > coex2/$filename.b.tran.sorted.txt
            join -t $'\t' -1 1 -2 1 coex2/RCA.b.tran.sorted.txt coex2/$filename.b.tran.sorted.txt > temp && mv temp coex2/RCA.b.tran.txt
            rm -rf coex2/$filename.b.tran.txt coex2/RCA.b.tran.sorted.txt coex2/$filename.b.tran.sorted.txt
        fi
    elif [ $(echo "$per2 > 50" | bc) -eq 1 ]; then
        # 只显示第一个满足条件的per2值
        if [ "$per2_shown" = false ] && [ -n "$per2" ]; then
            cut -f 1,4 $file > coex2/RCA.a.tran.txt
            awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/RCA.a.tran.txt > temp && mv temp coex2/RCA.a.tran.txt
            # 对RCA.tran.txt进行排序（保留表头）
            (head -n 1 coex2/RCA.a.tran.txt; tail -n +2 coex2/RCA.a.tran.txt | sort -k1,1) > temp && mv temp coex2/RCA.a.tran.txt
            per2_shown=true
        else
            cut -f 1,4 $file > coex2/$filename.a.tran.txt
            awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.a.tran.txt > temp && mv temp coex2/$filename.a.tran.txt
            # 对两个文件进行排序（保留表头）然后join
            (head -n 1 coex2/RCA.a.tran.txt; tail -n +2 coex2/RCA.a.tran.txt | sort -k1,1) > coex2/RCA.a.tran.sorted.txt
            (head -n 1 coex2/$filename.a.tran.txt; tail -n +2 coex2/$filename.a.tran.txt | sort -k1,1) > coex2/$filename.a.tran.sorted.txt
            join -t $'\t' -1 1 -2 1 coex2/RCA.a.tran.sorted.txt coex2/$filename.a.tran.sorted.txt > temp && mv temp coex2/RCA.a.tran.txt
            rm -rf coex2/$filename.a.tran.txt coex2/RCA.a.tran.sorted.txt coex2/$filename.a.tran.sorted.txt
        fi
    fi
done

for file in $(ls salmon3/Normal/*/quant.genes.sf)
do
    filename=$(awk -F'/' '{print $3}' <<< $file)
    path=$(ls RCA3/Normal/*/*/$filename*txt)
    p=$(awk -v var="$filename" -F',' '($1 == var || $11 == var) { print $30 }' Normal.csv | uniq)
    q=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk -F'[()]' '{print $2}')
    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR==5 {print $4}' $path)

    # 计算第四列第二行的值
    value2=$(awk 'NR==2 {print $4}' $path)

    # 计算第四列第三行的值
    value3=$(awk 'NR==3 {print $4}' $path)

    # 计算第四列第四行的值
    value4=$(awk 'NR==4 {print $4}' $path)

    # 计算第四列第三行到第四行的值的总和
    sum2=$(awk 'NR>=3 && NR<=4 {sum += $4} END {print sum}' $path)

    # 调试输出
    echo "Processing file: $path"
    echo "Values: sum=$sum, value2=$value2, value3=$value3, value4=$value4"

    # 确保 value2, value3, value4 和 sum 不是空的
    if [ -n "$value2" ] && [ -n "$value3" ] && [ -n "$value4" ] && [ -n "$sum" ]; then
        # 确保所有变量转换为普通浮点数格式
        sum=$(echo $sum | awk '{printf "%.10f\n", $1}')
        sum2=$(echo $sum2 | awk '{printf "%.10f\n", $1}')
        value2=$(echo $value2 | awk '{printf "%.10f\n", $1}')
        value3=$(echo $value3 | awk '{printf "%.10f\n", $1}')
        value4=$(echo $value4 | awk '{printf "%.10f\n", $1}')

        # 使用 bc 计算
        persum2=$(echo "scale=4; ($sum2/$sum)*100" | bc)
        per2=$(echo "scale=4; ($value2/$sum)*100" | bc)
        per3=$(echo "scale=4; ($value3/$sum)*100" | bc)
        per4=$(echo "scale=4; ($value4/$sum)*100" | bc)

        # 添加前导零
        persum2=$(printf "%.4f" $persum2)
        per2=$(printf "%.4f" $per2)
        per3=$(printf "%.4f" $per3)
        per4=$(printf "%.4f" $per4)
    else
        echo "Warning: Missing values in $path"
        per2=0
        per3=0
        per4=0
    fi

    if [ $(echo "$per2 < 50" | bc) -eq 1 ]; then
        cut -f 1,4 $file > coex2/$filename.b.tran.txt
        awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.b.tran.txt > temp && mv temp coex2/$filename.b.tran.txt
        # 对两个文件进行排序（保留表头）然后join
        (head -n 1 coex2/RCA.b.tran.txt; tail -n +2 coex2/RCA.b.tran.txt | sort -k1,1) > coex2/RCA.b.tran.sorted.txt
        (head -n 1 coex2/$filename.b.tran.txt; tail -n +2 coex2/$filename.b.tran.txt | sort -k1,1) > coex2/$filename.b.tran.sorted.txt
        join -t $'\t' -1 1 -2 1 coex2/RCA.b.tran.sorted.txt coex2/$filename.b.tran.sorted.txt > temp && mv temp coex2/RCA.b.tran.txt
        rm -rf coex2/$filename.b.tran.txt coex2/RCA.b.tran.sorted.txt coex2/$filename.b.tran.sorted.txt
    elif [ $(echo "$per2 > 50" | bc) -eq 1 ]; then
        cut -f 1,4 $file > coex2/$filename.a.tran.txt
        awk -v new_col="$q" 'BEGIN{FS=OFS="\t"} NR==1{$NF=new_col} {print}' coex2/$filename.a.tran.txt > temp && mv temp coex2/$filename.a.tran.txt
        # 对两个文件进行排序（保留表头）然后join
        (head -n 1 coex2/RCA.a.tran.txt; tail -n +2 coex2/RCA.a.tran.txt | sort -k1,1) > coex2/RCA.a.tran.sorted.txt
        (head -n 1 coex2/$filename.a.tran.txt; tail -n +2 coex2/$filename.a.tran.txt | sort -k1,1) > coex2/$filename.a.tran.sorted.txt
        join -t $'\t' -1 1 -2 1 coex2/RCA.a.tran.sorted.txt coex2/$filename.a.tran.sorted.txt > temp && mv temp coex2/RCA.a.tran.txt
        rm -rf coex2/$filename.a.tran.txt coex2/RCA.a.tran.sorted.txt coex2/$filename.a.tran.sorted.txt
    fi
done

awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex2/RCA.a.tran.txt | wc -l
awk 'NR==1 {for(i=1; i<=NF; i++) print $i}' coex2/RCA.b.tran.txt | wc -l # 查看文件的列数
