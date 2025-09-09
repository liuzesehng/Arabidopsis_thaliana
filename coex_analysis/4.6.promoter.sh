#!/bin/bash

# 使用samtools提取启动子序列
# 位置：NC_003071.7:16568745-16570745

# 设置变量
GENOME_FILE="/datapool/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_genomic.fa"
CHROM="NC_003071.7"
START="16568746"
END="16570745"
OUTPUT_PREFIX="rca_promoter.fasta"

# 检查基因组文件是否存在
if [ ! -f "$GENOME_FILE" ]; then
    echo "错误：基因组文件不存在 $GENOME_FILE"
    exit 1
fi

# 检查是否有samtools索引文件，如果没有则创建
if [ ! -f "${GENOME_FILE}.fai" ]; then
    echo "创建基因组索引文件..."
    samtools faidx "$GENOME_FILE"
fi

# 使用samtools faidx提取序列
echo "提取启动子序列：${CHROM}:${START}-${END}"
samtools faidx "$GENOME_FILE" "${CHROM}:${START}-${END}" > "${OUTPUT_PREFIX}"

# 检查提取是否成功
if [ $? -eq 0 ]; then
    echo "启动子序列成功提取到：${OUTPUT_PREFIX}"
    echo "序列信息："
    head -n 2 "${OUTPUT_PREFIX}"
    echo "序列长度：$(tail -n +2 "${OUTPUT_PREFIX}" | tr -d '\n' | wc -c) bp"
else
    echo "错误：序列提取失败"
    exit 1
fi

# 显示输出文件的完整路径
echo "输出文件完整路径：$(pwd)/${OUTPUT_PREFIX}"
