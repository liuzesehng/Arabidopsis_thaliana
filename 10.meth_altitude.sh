#!/bin/bash



# 定义处理每个文件对或单文件的函数
process_bed() {

    local files="$1"
    local file="$2"
    filename=${files/.txt/}
    
    # 提取数据
    while IFS=$'\t' read -r col1 start end; do


        # 定义处理单个文件的函数
        process_file() {
            local i="$1"
            local col1="$2"
            local start="$3"
            local end="$4"
            local file="$5"
            local filename="$6"
            # 使用 awk 打印符合条件的第四列的值
            awk -v col1="$col1" -v start="$start" -v end="$end" '$1 == col1 && $2 == start && $3 == end {print $4}' "$i" >> meth/$file.$filename.me_count.txt           
        }

        # 导出函数以便 GNU Parallel 使用
        export -f process_file

        # 并行处理每个文件，并传递所需参数
        grep "$file" meth/"$files" | parallel -j 90 process_file {} "$col1" "$start" "$end" "$file" "$filename"
        # 计算均值
        mean=$(awk '{sum+=$1} END {print sum/NR}' meth/$file.$filename.me_count.txt) 
        # 计算标准差
        std_dev=$(awk -v mean=$mean '{sum+=($1-mean)^2} END {print sqrt(sum/NR)}' meth/$file.$filename.me_count.txt)
        # 计算变异系数
        cv=$(awk -v std_dev=$std_dev -v mean=$mean 'BEGIN {print (std_dev/mean) * 100}')

        if [ -z "$cv" ]; then
            cv=0
        fi

        echo -e "NC_003071.7\t$start\t$end\t$cv" >> meth/$file.$filename.CV.bedGraph

        rm -rf meth/$file.$filename.me_count.txt


    done < meth/RCA.up_down.bed
}

# 导出函数以便 GNU Parallel 使用
export -f process_bed

# 使用 parallel 处理所有文件类型
#parallel -j 6 process_bed ::: low_altitude.txt mid_altitude.txt ::: CG CHG CHH
parallel -j 1 process_bed ::: low_altitude.txt ::: CHH
exit 0

for file in $(awk -F'\t' '$4 != "" && $4 > 0 && $4 <= 1000 {print $1}' Ala.ab.location.txt)
do
    altitude=$(awk -F'\t' -v var="$file" '$1 == var && $2 != "" {print $2}' meth/Ala.abnor.meth.txt)
    for i in $(find Abnormal_me/bismark -name "$altitude".*deduplicated*.bedGraph | grep -v "RCA") 
    do
        echo "$i" >> meth/low_altitude.txt
    done  
done

for file in $(awk -F'\t' '$4 != "" && $4 > 0 && $4 <= 1000 {print $1}' Ala.normal.location.txt)
do
    altitude=$(awk -F'\t' -v var="$file" '$1 == var && $2 != "" {print $2}' meth/Ala.nor.meth.txt)
    for i in $(find Normal_me/bismark -name "$altitude".*deduplicated*.bedGraph | grep -v "RCA")
    do
        echo "$i" >> meth/low_altitude.txt
    done
done

for file in $(awk -F'\t' '$4 != "" && $4 > 1000 && $4 <= 3500 {print $1}' Ala.ab.location.txt)
do
    altitude=$(awk -F'\t' -v var="$file" '$1 == var && $2 != "" {print $2}' meth/Ala.abnor.meth.txt)
    for i in $(find Abnormal_me/bismark -name "$altitude".*deduplicated*.bedGraph | grep -v "RCA")
    do
        echo "$i" >> meth/mid_altitude.txt
    done 
done

for file in $(awk -F'\t' '$4 != "" && $4 > 1000 && $4 <= 3500 {print $1}' Ala.normal.location.txt)
do
    altitude=$(awk -F'\t' -v var="$file" '$1 == var && $2 != "" {print $2}' meth/Ala.nor.meth.txt)
    for i in $(find Normal_me/bismark -name "$altitude".*deduplicated*.bedGraph | grep -v "RCA")
    do
        echo "$i" >> meth/mid_altitude.txt
    done  
done