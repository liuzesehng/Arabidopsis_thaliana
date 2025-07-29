#!/bin/bash

# 定义处理每个文件对或单文件的函数
process_bed() {
    file="$1"

    # 提取数据
    while IFS=$'\t' read -r col1 start end; do

        # 定义处理单个文件的函数
            process_file() {
            local i="$1"
            local col1="$2"
            local start="$3"
            local end="$4"
            local file="$5"
            # 使用 awk 打印符合条件的第四列的值
            awk -v col1="$col1" -v start="$start" -v end="$end" '$1 == col1 && $2 == start && $3 == end {print $4}' "$i" >> meth/$file.me_count.txt           
        }

        # 导出函数以便 GNU Parallel 使用
        export -f process_file

        # 并行处理每个文件，并传递所需参数
        grep "$file" meth/meth.txt | parallel -j 84 process_file {} "$col1" "$start" "$end" "$file"


        # 计算均值
        mean=$(awk '{sum+=$1} END {print sum/NR}' meth/$file.me_count.txt)
        # 计算标准差
        std_dev=$(awk -v mean=$mean '{sum+=($1-mean)^2} END {print sqrt(sum/NR)}' meth/$file.me_count.txt)
        # 计算变异系数
        cv=$(awk -v std_dev=$std_dev -v mean=$mean 'BEGIN {print (std_dev/mean) * 100}')

        if [ -z "$mean" ];then
            mean=0
        fi

        if [ -z "$cv" ];then
            cv=0
        fi

        echo -e "NC_003071.7\t$start\t$end\t$mean\t$cv" >> meth/$file.CV.bedGraph

        rm -rf meth/$file.me_count.txt

    done < meth/RCA.up_down.bed
}

# 导出函数以便 GNU Parallel 使用
export -f process_bed

# 使用 parallel 处理所有文件类型
parallel -j 3 process_bed ::: CG CHG CHH

exit 0 

# 提取数据
for file in $(find Abnormal_me/bismark Normal_me/bismark -name *RCA.promoter_terminator.processed.bedGraph)
do
    echo "$file" >> meth/meth.txt
done

