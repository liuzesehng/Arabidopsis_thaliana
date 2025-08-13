#!/bin/bash

# 输入和输出文件路径
INPUT_FILE="total/CytoscapeInput-edges-rca.csv"
OUTPUT_FILE="total/gene_id.tsv"

# 检查输入文件是否存在
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误：找不到输入文件 $INPUT_FILE"
    exit 1
fi

# 处理CSV文件
echo "正在处理文件..."

# 使用awk处理文件
awk -F',' '
BEGIN {
    OFS="\t"
    print "name", "shared_name", "weight"
}
NR > 1 {
    # 去除引号
    gsub(/"/, "", $2)
    gsub(/"/, "", $6)
    
    name = $2
    weight = $6
    
    # 分割 name 列
    if (index(name, " (interacts with) ") > 0) {
        split(name, parts, " \\(interacts with\\) ")
        gene1 = parts[1]
        gene2 = parts[2]
        
        # AT2G39730 放在第一列
        if (gene1 == "AT2G39730") {
            first_col = gene1
            second_col = gene2
        } else if (gene2 == "AT2G39730") {
            first_col = "AT2G39730"
            second_col = gene1
        } else {
            first_col = gene1
            second_col = gene2
        }
        
        print first_col, second_col, weight
    }
}' "$INPUT_FILE" | sort -k3,3nr > "$OUTPUT_FILE"

echo "处理完成！输出文件：$OUTPUT_FILE"
echo "共处理 $(tail -n +2 "$OUTPUT_FILE" | wc -l) 行数据"
