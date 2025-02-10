#!/bin/bash
#bismark_genome_preparation --parallel 16 /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen

# 定义处理每个文件对或单文件的函数
process_fastq() {
    local p1="$1"
    
    if [[ $p1 =~ ^Abnormal_me/.*_1\.fit\.fastq\.gz$ || $p1 =~ ^Normal_me/.*_1\.fit\.fastq\.gz$ ]]; then
		pref=${p1/_1.fit.fastq.gz/}
		ref=$(dirname $(dirname "$p1"))
		ref2=$(basename "$pref")
		p2=$pref"_2.fit.fastq.gz"

        if [ -d $ref/bismark/$ref2 ]; then
            echo "dir exists"
        else
            mkdir -p $ref/bismark/$ref2
        fi
		bismark /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/bismark -p 16 -o $ref/bismark/$ref2 -1 $p1 -2 $p2
		samtools sort -n -@ 16 -o $ref/bismark/$ref2/$ref2.nameSorted.bam $ref/bismark/$ref2/$ref2"_1.fit_bismark_bt2_pe.bam"
		deduplicate_bismark -p --bam --output_dir $ref/bismark/$ref2 $ref/bismark/$ref2/$ref2.nameSorted.bam

    else
		pref=${p1/.fit.fastq.gz/}
		ref=$(dirname $(dirname "$p1"))
		ref2=$(basename "$pref")

        if [ -d $ref/bismark/$ref2 ]; then
            echo "dir exists"
        else
            mkdir -p $ref/bismark/$ref2
        fi
		bismark /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/bismark -p 16 -o $ref/bismark/$ref2 $p1
		samtools sort -n -@ 16 -o $ref/bismark/$ref2/$ref2.nameSorted.bam $ref/bismark/$ref2/$ref2".fit_bismark_bt2.bam"
		deduplicate_bismark -s --bam --output_dir $ref/bismark/$ref2 $ref/bismark/$ref2/$ref2.nameSorted.bam
    fi
}

# 导出函数以便 GNU Parallel 使用
export -f process_fastq

# 找出所有需要处理的文件

parallel -j 4 process_fastq ::: $(find Abnormal_me Normal_me -name "*_1.fit.fastq.gz" | sed -n '1,300p')	
#parallel -j 4 process_fastq ::: $(find Normal_me -name "*.fit.fastq.gz" | grep -v -E "_1|_2")
