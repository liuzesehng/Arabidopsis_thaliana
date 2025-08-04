#!/bin/bash
# 输入和输出文件路径
INPUT_FILE="total/CytoscapeInput-edges-turquoise.degree.node.csv"
OUTPUT_FILE="total/gene_id.tsv"

# 检查输入文件是否存在
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误：找不到输入文件 $INPUT_FILE"
    exit 1
fi

# 处理CSV文件
echo "正在处理文件..."

# 使用awk处理文件
# 输出表头到文件
awk -F',' 'BEGIN { OFS="\t" } NR == 1 { print $10, $5 }' "$INPUT_FILE" > "$OUTPUT_FILE"

# 追加排序后的数据
awk -F',' 'BEGIN { OFS="\t" } NR > 1 && NR < 181 { print $10, $5 }' "$INPUT_FILE" | \
sort -k2,2nr >> "$OUTPUT_FILE"

echo "处理完成！输出文件：$OUTPUT_FILE"
echo "共处理 $(tail -n +2 "$OUTPUT_FILE" | wc -l) 行数据"