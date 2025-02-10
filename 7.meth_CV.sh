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

        if [ -z "$cv" ];then
            cv=0
        fi

        echo -e "NC_003071.7\t$start\t$end\t$cv" >> meth/$file.$filename.CV.bedGraph

        rm -rf meth/$file.$filename.me_count.txt

    done < meth/RCA.up_down.bed
}

# 导出函数以便 GNU Parallel 使用
export -f process_bed

# 使用 parallel 处理所有文件类型
#parallel -j 9 process_bed ::: 10C.txt 16C.txt 22C.txt ::: CG CHG CHH
parallel -j 1 process_bed ::: 22C.txt ::: CHH

exit 0

for file in $(find Abnormal_me/bismark -name *.nameSorted.deduplicated*.bedGraph | grep -v "RCA")
do
    ref=$(awk -F'/' '{print $3}' <<< $file)
    ref2=$(awk -v var="$ref" -F',' '($1 == var || $11 == var) { print $30 }' Abnormal_me/Abnormal.me.csv | uniq) 
    tem=$(cat meth/biosample_result_me.txt | grep -A 9 "$ref2" | grep 'growth temperature="' | sed 's/.*\/growth temperature="\([^"]*\)".*/\1/')

    if [ $tem == "10C" ]; then
        echo "$file" >> meth/10C.txt
    else
        echo "$file" >> meth/16C.txt
    fi
done

for file in $(find Normal_me/bismark -name *.nameSorted.deduplicated*.bedGraph | grep -v "RCA")
do
    echo "$file" >> meth/22C.txt
done

# 输出文件名
output_file="meth/RCA.up_down.bed"

# 起始值和结束值
start=16568745
end=16575692

# 打开输出文件
exec > "$output_file"

# 生成内容
for ((i=start; i<end; i++)); do
    echo -e "NC_003071.7\t$i\t$((i+1))"
done

