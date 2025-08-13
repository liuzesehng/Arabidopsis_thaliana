#!/bin/bash

# 定义处理每个文件对或单文件的函数
process_bam() {
    local bam="$1"
    
	pref=${bam/.bam/}
	ref=$(dirname "$bam")
	ref2=$(basename "$pref")

	p1=$(ls $ref/*E_report.txt)
	if [[ $p1 =~ ^Abnormal_me/.*fit_bismark_bt2_PE_report\.txt$ || $p1 =~ ^Normal_me/.*fit_bismark_bt2_PE_report\.txt$ ]]; then

		bismark_methylation_extractor -p --no_overlap --comprehensive --cytosine_report --gzip --counts --bedGraph --genome_folder /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/bismark -o $ref --buffer_size 60G --parallel 10 --CX $bam

	else
		bismark_methylation_extractor -s --no_overlap --comprehensive --cytosine_report --gzip --counts --bedGraph --genome_folder /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/bismark -o $ref --buffer_size 60G --parallel 10 --CX $bam

	fi

	zcat $pref.CX_report.txt.gz | awk '$4+$5>=3' | grep -P "\tCG\t" | gzip - > $pref.CG.CX_report.txt.gz
	zcat $pref.CX_report.txt.gz | awk '$4+$5>=3' | grep -P "\tCHG\t" | gzip - > $pref.CHG.CX_report.txt.gz
	zcat $pref.CX_report.txt.gz | awk '$4+$5>=3' | grep -P "\tCHH\t" | gzip - > $pref.CHH.CX_report.txt.gz
		
	bismark2bedGraph --CX --buffer_size 60G --dir $ref -o $ref2.CG $ref/"CpG_context_"$ref2.txt.gz
	bismark2bedGraph --CX --buffer_size 60G --dir $ref -o $ref2.CHG $ref/"CHG_context_"$ref2.txt.gz
	bismark2bedGraph --CX --buffer_size 60G --dir $ref -o $ref2.CHH $ref/"CHH_context_"$ref2.txt.gz
		
	zcat $pref.CG.gz > $pref.CG.bedgraph
	zcat $pref.CHG.gz > $pref.CHG.bedgraph
	zcat $pref.CHH.gz > $pref.CHH.bedgraph
		
	bedGraphToBigWig $pref.CG.bedgraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt $pref.CG.bw
	bedGraphToBigWig $pref.CHG.bedgraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt $pref.CHG.bw
	bedGraphToBigWig $pref.CHH.bedgraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt $pref.CHH.bw

	rm -rf $ref/"CpG_context_"$ref2.txt.gz $ref/"CHG_context_"$ref2.txt.gz $ref/"CHH_context_"$ref2.txt.gz $pref.CG.gz $pref.CHG.gz $pref.CHH.gz $pref.CG.gz.bismark.cov.gz $pref.CHG.gz.bismark.cov.gz $pref.CHH.gz.bismark.cov.gz $pref.CG.bedgraph $pref.CHG.bedgraph $pref.CHH.bedgraph $pref.bedGraph.gz $pref.bismark.cov.gz
}

# 导出函数以便 GNU Parallel 使用
export -f process_bam

# 找出所有需要处理的文件

parallel -j 2 process_bam ::: $(find Abnormal_me/bismark Normal_me/bismark -name *deduplicated.bam | sed -n '201,400p')

