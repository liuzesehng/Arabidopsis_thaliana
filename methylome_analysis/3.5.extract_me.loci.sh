#!/bin/bash
echo -e "NC_003071.7\t16570745\t16573692" > gene_RCA
echo -e "NC_003071.7\t16568745\t16570745" >> gene_RCA_promoter
echo -e "NC_003071.7\t16573692\t16575692" >> gene_RCA_terminator
# 定义处理每个文件对或单文件的函数
process_bw() {
    local bw="$1"
    
	pref=${bw/.bw/}

    bigWigToBedGraph $bw $pref.bedGraph &&
    #awk '$1 == "NC_003071.7" && $2 >= 16568745 && $2 <= 16570745 && $3 <= 16570745' $pref.bedGraph > $pref.RCA.promoter.bedGraph
    #awk '$1 == "NC_003071.7" && $2 >= 16570745 && $2 <= 16573692 && $3 <= 16573692' $pref.bedGraph > $pref.RCA.bedGraph
    #awk '$1 == "NC_003071.7" && $2 >= 16573692 && $2 <= 16575692 && $3 <= 16575692' $pref.bedGraph > $pref.RCA.terminator.bedGraph
    awk '$1 == "NC_003071.7" && $2 >= 16568745 && $2 <= 16575692 && $3 <= 16575692' $pref.bedGraph > $pref.RCA.promoter_terminator.bedGraph

}

# 导出函数以便 GNU Parallel 使用
export -f process_bw

# 找出所有需要处理的文件

parallel -j 64 process_bw ::: $(find Abnormal_me/bismark Normal_me/bismark -name *.bw)

# 定义处理每个文件对或单文件的函数
process_bed() {
    local file="$1"
    
	pref=${file/.bedGraph/}

    # 读取输入文件的每一行
    while IFS=$'\t' read -r chr start end value; do
    # 计算区间长度
    length=$((end - start))
    
    if [ "$length" -gt 1 ]; then
        
        # 拆分区间
        for ((i=0; i<length; i++)); do
        new_start=$((start + i))
        new_end=$((new_start + 1))
        if [[ "$value" == *.* ]]; then
            echo -e "$chr\t$new_start\t$new_end\t$value" | awk '{printf "%s\t%d\t%d\t%.4f\n", $1, $2, $3, $4}' >> "$pref.processed.bedGraph"
        else
            echo -e "$chr\t$new_start\t$new_end\t$value" | awk '{printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4}' >> "$pref.processed.bedGraph"
        fi
        done
    else
        # 如果区间长度等于1，直接写入输出文件
        if [[ "$value" == *.* ]]; then
            echo -e "$chr\t$start\t$end\t$value" | awk '{printf "%s\t%d\t%d\t%.4f\n", $1, $2, $3, $4}' >> "$pref.processed.bedGraph"
        else
            echo -e "$chr\t$start\t$end\t$value" | awk '{printf "%s\t%d\t%d\t%d\n", $1, $2, $3, $4}' >> "$pref.processed.bedGraph"
        fi
    fi
    done < "$file"


}

# 导出函数以便 GNU Parallel 使用
export -f process_bed

# 找出所有需要处理的文件

parallel -j 90 process_bed ::: $(find Abnormal_me/bismark Normal_me/bismark -name *RCA.promoter_terminator.bedGraph)
exit 0

# 提取数据
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

