#!/bin/bash


# 定义处理每个文件对或单文件的函数
process_fastq() {
    local p1="$1"
    
    pref=${p1/R1.fastq.gz/}
    ref=$(basename $pref)
    p2=$pref"R2.fastq.gz"

    dir=$(echo $p1 | cut -d '/' -f 2)

    fit1=methdata/$dir/$dir".R1.fit.fastq.gz"
    fit2=methdata/$dir/$dir".R2.fit.fastq.gz"
    fastp -z 9 -w 4 -i $p1 -I $p2 -o $fit1 -O $fit2 -j methdata/$dir/$dir"_fastp.json" -h methdata/$dir/$dir"_fastp.html"
}

# 导出函数以便 GNU Parallel 使用
export -f process_fastq

# 找出所有需要处理的文件

parallel -j 16 process_fastq ::: $(find methdata -name "*R1.fastq.gz")
