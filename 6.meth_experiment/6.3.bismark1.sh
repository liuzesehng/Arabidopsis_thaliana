#!/bin/bash
# 定义处理每个文件对或单文件的函数
process_fastq() {
    local p1="$1"

		pref=${p1/R1.fit.fastq.gz/}
		ref=$(echo $pref | cut -d '/' -f 2)
		p2=$pref"R2.fit.fastq.gz"

        if [ -d methdata/bismark/$ref ]; then
            echo "dir exists"
        else
            mkdir -p methdata/bismark/$ref
        fi

		bismark /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/bismark -p 16 -o methdata/bismark/$ref -1 $p1 -2 $p2
		samtools sort -n -@ 16 -o methdata/bismark/$ref/$ref.nameSorted.bam methdata/bismark/$ref/$ref".R1.fit_bismark_bt2_pe.bam"
		deduplicate_bismark -p --bam --output_dir methdata/bismark/$ref methdata/bismark/$ref/$ref.nameSorted.bam
}

# 导出函数以便 GNU Parallel 使用
export -f process_fastq

# 找出所有需要处理的文件

parallel -j 4 process_fastq ::: $(find methdata -name "*R1.fit.fastq.gz" | sed -n '1,12p')	
