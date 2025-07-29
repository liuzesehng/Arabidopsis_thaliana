#!/bin/bash
# 定义处理每个文件对或单文件的函数
process_bed() {
    local file="$1"

    pref=$(basename $file)
    pref=${pref/.CX_report.txt.gz}
    pref=${pref/.nameSorted.deduplicated}
    python /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/program/cgmaptools-0.1.3/bin/BismarkToCGmap -i <(zcat $file) -o $pref.CGmap
    python /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/program/cgmaptools-0.1.3/bin/CGmapToRegion -i $pref.CGmap -r gene_RCA > $pref.rca.mtr
    python /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/program/cgmaptools-0.1.3/bin/CGmapToRegion -i $pref.CGmap -r gene_RCA_promoter > $pref.rca_promoter.mtr
    python /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/program/cgmaptools-0.1.3/bin/CGmapToRegion -i $pref.CGmap -r gene_RCA_terminator > $pref.rca_terminator.mtr
    count1=$(cut -f 6 $pref.rca.mtr)
    count2=$(cut -f 6 $pref.rca_promoter.mtr)
    count3=$(cut -f 6 $pref.rca_terminator.mtr)
    echo -e "$pref.rca\t$count1" >> meth/methylation_count.txt
    echo -e "$pref.rca_promoter\t$count2" >> meth/methylation_count.txt
    echo -e "$pref.rca_terminator\t$count3" >> meth/methylation_count.txt
    rm -rf $pref.CGmap $pref.rca.mtr $pref.rca_promoter.mtr $pref.rca_terminator.mtr
}

# 导出函数以便 GNU Parallel 使用
export -f process_bed

# 使用 parallel 处理所有文件类型
parallel -j 64 process_bed ::: $(ls */bismark/*/*report.txt.gz | grep -v ".deduplicated.CX")

exit 0

for file in $(ls */bismark/*/*RCA*bedGraph)
do
    file_name=$(basename $file)
    file_name=${file_name/.bedGraph} 
    file_name=${file_name/.nameSorted.deduplicated}

    #count=$(awk '{sum += ($3 - $2) * $4} END {print sum}' $file)
    count=$(awk '{sum += ($3 - $2) * $4; diff_sum += ($3 - $2)} END {if (diff_sum > 0) print sum / diff_sum; else print "Division by zero error"}' $file)

    echo -e "$file_name\t$count" >> meth/methylation_count.txt
done
#find RCA/Abnormal RCA/Normal -type f -name "*.txt"



