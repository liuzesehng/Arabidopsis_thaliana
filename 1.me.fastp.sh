#!/bin/bash


# 定义处理每个文件对或单文件的函数
process_fastq() {
    local p1="$1"
    
    if [[ $p1 =~ ^Abnormal_me/.*_1\.fastq\.gz$ || $p1 =~ ^Normal_me/.*_1\.fastq\.gz$ ]]; then
        pref=${p1/_1.fastq.gz/}
        ref=$(basename $pref)
        p2=$pref"_2.fastq.gz"

        if [ -d $pref ]; then
            echo "dir exists"
        else
            mkdir -p $pref
        fi

        fit1=$pref/$ref"_1.fit.fastq.gz"
        fit2=$pref/$ref"_2.fit.fastq.gz"
        fastp -z 9 -w 4 -i $p1 -I $p2 -o $fit1 -O $fit2 -j $pref/$ref"_fastp.json" -h $pref/$ref"_fastp.html"
    else
        pref=${p1/.fastq.gz/}
        ref=$(basename $pref)

        if [ -d $pref ]; then
            echo "dir exists"
        else
            mkdir -p $pref
        fi

        fit=$pref/$ref".fit.fastq.gz"
        fastp -z 9 -w 4 -i $p1 -o $fit -j $pref/$ref"_fastp.json" -h $pref/$ref"_fastp.html"
    fi
}

# 导出函数以便 GNU Parallel 使用
export -f process_fastq

# 找出所有需要处理的文件

parallel -j 16 process_fastq ::: $(find Abnormal_me Normal_me -name "*_1.fastq.gz")

parallel -j 16 process_fastq ::: $(find Normal_me -name "*.fastq.gz" | grep -v -E "_1|_2")

