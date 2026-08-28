#!/bin/bash
# # 定义处理每个文件对或单文件的函数
# process_cx() {
#     local file="$1"
    
# 	pref=${file/.txt/}
#     awk '$1 == "NC_003071.7" && $2 >= 16568746 && $2 <= 16575692' $file > $pref.RCA.promoter_terminator.txt

# }

# # 导出函数以便 GNU Parallel 使用
# export -f process_cx

# # 找出所有需要处理的文件
# parallel -j 64 process_cx ::: $(ls */bismark/*/*CX_report.txt)

# 初始化列名数组
columns=("Temperature")

# 动态生成 CG_1 到 CG_262
for i in {1..262}; do
    columns+=("CG_$i")
done

# 动态生成 CHG_1 到 CHG_337
for i in {1..337}; do
    columns+=("CHG_$i")
done

# 动态生成 CHH_1 到 CHH_1826
for i in {1..1826}; do
    columns+=("CHH_$i")
done

# 动态生成 promoter、terminator数据
# columns+=("CG_promoter" "CG" "CG_terminator" "CHG_promoter" "CHG" "CHG_terminator" "CHH_promoter" "CHH" "CHH_terminator")

# 动态生成 snp_1 到 snp_104
for i in {1..104}; do
    columns+=("snp_$i")
done

columns+=("total" "α" "β" "β1" "β2" "α/%" "β/%" "β1/%" "β2/%")
# 写入表头
echo -e "$(IFS=$'\t'; echo "${columns[*]}")" > RCA/Alt.weighted.snp_meth.tsv

# 定义函数处理 cx 文件
process_cx_coverage() {
    local prefix=$1   # CG / CHG / CHH
    local max=$2
    local folder=$3   # *.CX_report.txt.gz

    for i in $(seq 1 $max); do
        j=$(awk "NR == $i {print \$3}" meth/${prefix}.CV.before.bedGraph)

        value=$(awk -v pos="$j" '$2 == pos {print $4 + $5; exit}' "$folder")
        if [ -z "$value" ]; then
            value=""
        fi
        col+=("$value")
    done
}

# 定义函数处理 cx1 文件
process_cx1_coverage() {
    local prefix=$1  # CG, CHG, CHH
    local max=$2     # 最大迭代数
    local folder=$3  # 文件夹路径
    for i in $(seq 1 $max); do
        j=$(awk "NR == $i {print \$3}" meth/${prefix}.CV.before.bedGraph)
        value=$(awk -v pos="$j" '$2 == pos {print $4 + $5; exit}' "$folder")
        if [ -z "$value" ]; then
            value=""
        fi
        col1+=("$value")
    done
}

# 定义函数处理 cx2 文件
process_cx2_coverage() {
    local prefix=$1  # CG, CHG, CHH
    local max=$2     # 最大迭代数
    local folder=$3  # 文件夹路径
    for i in $(seq 1 $max); do
        j=$(awk "NR == $i {print \$3}" meth/${prefix}.CV.before.bedGraph)
        value=$(awk -v pos="$j" '$2 == pos {print $4 + $5; exit}' "$folder")
        if [ -z "$value" ]; then
            value=""
        fi
        col2+=("$value")
    done
}
#RCA数据整理
for path in $(ls RCA/Abnormal/*/*.txt)
do    
    # 使用cut命令和斜杠（/）作为分隔符，提取第3个字段
    tem=$(echo $path | cut -d'/' -f3)
    tem=${tem/C}

    # 提取名称
    me=${path/.RCA.transcript.txt/}
    name=$(basename ${me})

    col=()
    col1=()
    col2=()

    name2=$(awk '$1 == "'$name'" {print $2}' meth/Ala.abnor.meth.txt)

    # 将 name2 拆分成数组
    name_array=($name2)

    # 判断是否有两个或更多的值
    if [ ${#name_array[@]} -ge 2 ]; then
        col1+=("$tem")
        col2+=("$tem")

        # 处理 CG, CHG, CHH 数据
        process_cx1_coverage "CG" 262 "Abnormal_me/bismark/${name_array[0]}/${name_array[0]}.nameSorted.deduplicated.CG.CX_report.RCA.promoter_terminator.txt"
        process_cx1_coverage "CHG" 337 "Abnormal_me/bismark/${name_array[0]}/${name_array[0]}.nameSorted.deduplicated.CHG.CX_report.RCA.promoter_terminator.txt"
        process_cx1_coverage "CHH" 1826 "Abnormal_me/bismark/${name_array[0]}/${name_array[0]}.nameSorted.deduplicated.CHH.CX_report.RCA.promoter_terminator.txt"
        process_cx2_coverage "CG" 262 "Abnormal_me/bismark/${name_array[1]}/${name_array[1]}.nameSorted.deduplicated.CG.CX_report.RCA.promoter_terminator.txt"
        process_cx2_coverage "CHG" 337 "Abnormal_me/bismark/${name_array[1]}/${name_array[1]}.nameSorted.deduplicated.CHG.CX_report.RCA.promoter_terminator.txt"
        process_cx2_coverage "CHH" 1826 "Abnormal_me/bismark/${name_array[1]}/${name_array[1]}.nameSorted.deduplicated.CHH.CX_report.RCA.promoter_terminator.txt"

        # 处理 snp 数据
        p=$(awk -v var="$name" -F',' '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)

        q=$(cat biosample_result.txt | grep -A 9 $p | grep "accession number=" | sed 's/.*\/accession number="\([^"]*\)".*/\1/')

        # snp_1 到 snp_104
        for snp in {1..104}
        do
            j=$(awk -v i="$snp" 'NR > 1 && NR == i + 1 {print $2}' snp/filtered_snp_matrix.csv)
            eval snp_$snp=$(awk -v q="$q" -v j="$j" 'NR==1 {for(i=1;i<=NF;i++) if($i==q) {col=i; break}; next} col && $2==j {print $col}' snp/filtered_snp_matrix.csv)
            if [ -z "$(eval echo '$snp_'$snp)" ]; then
                col1+=("")
                col2+=("")
            else
                col1+=($(eval echo '$snp_'$snp))
                col2+=($(eval echo '$snp_'$snp))
            fi
        done

        # 计算第四列第二行到第四行的值的总和
        # sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)
        sum=$(awk 'NR==5 {print $4}' $path)

        # 计算第四列第二行的值
        value2=$(awk 'NR==2 {print $4}' $path)

        # 计算第四列第三行的值
        value3=$(awk 'NR==3 {print $4}' $path)

        # 计算第四列第四行的值
        value4=$(awk 'NR==4 {print $4}' $path)

        # 计算第四列第三行到第四行的值的总和
        sum2=$(awk 'NR>=3 && NR<=4 {sum += $4} END {print sum}' $path)

        # 调试输出
        echo "Processing file: $path"
        echo "Values: sum=$sum, value2=$value2, value3=$value3, value4=$value4"

        # 确保 value2, value3, value4 和 sum 不是空的
        if [ -n "$value2" ] && [ -n "$value3" ] && [ -n "$value4" ] && [ -n "$sum" ]; then
            # 确保所有变量转换为普通浮点数格式
            sum=$(echo $sum | awk '{printf "%.10f\n", $1}')
            sum2=$(echo $sum2 | awk '{printf "%.10f\n", $1}')
            value2=$(echo $value2 | awk '{printf "%.10f\n", $1}')
            value3=$(echo $value3 | awk '{printf "%.10f\n", $1}')
            value4=$(echo $value4 | awk '{printf "%.10f\n", $1}')

            # 使用 bc 计算
            persum2=$(echo "scale=4; ($sum2/$sum)*100" | bc)
            per2=$(echo "scale=4; ($value2/$sum)*100" | bc)
            per3=$(echo "scale=4; ($value3/$sum)*100" | bc)
            per4=$(echo "scale=4; ($value4/$sum)*100" | bc)

            # 添加前导零
            persum2=$(printf "%.4f" $persum2)
            per2=$(printf "%.4f" $per2)
            per3=$(printf "%.4f" $per3)
            per4=$(printf "%.4f" $per4)
        else
            echo "Warning: Missing values in $path"
            per2=0
            per3=0
            per4=0
        fi
        col1+=($sum $value2 $sum2 $value3 $value4 $per2 $persum2 $per3 $per4)
        col2+=($sum $value2 $sum2 $value3 $value4 $per2 $persum2 $per3 $per4)

        echo -e "$(IFS=$'\t'; echo "${col1[*]}")" >> RCA/Alt.snp_meth.tsv
        echo -e "$(IFS=$'\t'; echo "${col2[*]}")" >> RCA/Alt.snp_meth.tsv
        col1=()
        col2=()

    else
        col+=("$tem")
        if [ -z "$name2" ]; then
            for i in {1..2425}
            do
                col+=("")
            done
        else
            # 处理 CG, CHG, CHH 数据
            process_cx_coverage "CG" 262 "Abnormal_me/bismark/$name2/$name2.nameSorted.deduplicated.CG.CX_report.RCA.promoter_terminator.txt"
            process_cx_coverage "CHG" 337 "Abnormal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHG.CX_report.RCA.promoter_terminator.txt"
            process_cx_coverage "CHH" 1826 "Abnormal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHH.CX_report.RCA.promoter_terminator.txt"
        fi

        # 处理 snp 数据
        p=$(awk -v var="$name" -F',' '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)

        q=$(cat biosample_result.txt | grep -A 9 $p | grep "accession number=" | sed 's/.*\/accession number="\([^"]*\)".*/\1/')

        # snp_1 到 snp_104
        for snp in {1..104}
        do
            j=$(awk -v i="$snp" 'NR > 1 && NR == i + 1 {print $2}' snp/filtered_snp_matrix.csv)
            eval snp_$snp=$(awk -v q="$q" -v j="$j" 'NR==1 {for(i=1;i<=NF;i++) if($i==q) {col=i; break}; next} col && $2==j {print $col}' snp/filtered_snp_matrix.csv)
            if [ -z "$(eval echo '$snp_'$snp)" ]; then
                col+=("")
            else
                col+=($(eval echo '$snp_'$snp))
            fi
        done

        # 计算第四列第二行到第四行的值的总和
        # sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)
        sum=$(awk 'NR==5 {print $4}' $path)

        # 计算第四列第二行的值
        value2=$(awk 'NR==2 {print $4}' $path)

        # 计算第四列第三行的值
        value3=$(awk 'NR==3 {print $4}' $path)

        # 计算第四列第四行的值
        value4=$(awk 'NR==4 {print $4}' $path)

        # 计算第四列第三行到第四行的值的总和
        sum2=$(awk 'NR>=3 && NR<=4 {sum += $4} END {print sum}' $path)

        # 调试输出
        echo "Processing file: $path"
        echo "Values: sum=$sum, value2=$value2, value3=$value3, value4=$value4"

        # 确保 value2, value3, value4 和 sum 不是空的
        if [ -n "$value2" ] && [ -n "$value3" ] && [ -n "$value4" ] && [ -n "$sum" ]; then
            # 确保所有变量转换为普通浮点数格式
            sum=$(echo $sum | awk '{printf "%.10f\n", $1}')
            sum2=$(echo $sum2 | awk '{printf "%.10f\n", $1}')
            value2=$(echo $value2 | awk '{printf "%.10f\n", $1}')
            value3=$(echo $value3 | awk '{printf "%.10f\n", $1}')
            value4=$(echo $value4 | awk '{printf "%.10f\n", $1}')

            # 使用 bc 计算
            persum2=$(echo "scale=4; ($sum2/$sum)*100" | bc)
            per2=$(echo "scale=4; ($value2/$sum)*100" | bc)
            per3=$(echo "scale=4; ($value3/$sum)*100" | bc)
            per4=$(echo "scale=4; ($value4/$sum)*100" | bc)

            # 添加前导零
            persum2=$(printf "%.4f" $persum2)
            per2=$(printf "%.4f" $per2)
            per3=$(printf "%.4f" $per3)
            per4=$(printf "%.4f" $per4)
        else
            echo "Warning: Missing values in $path"
            per2=0
            per3=0
            per4=0
        fi
        col+=($sum $value2 $sum2 $value3 $value4 $per2 $persum2 $per3 $per4)

        echo -e "$(IFS=$'\t'; echo "${col[*]}")" >> RCA/Alt.weighted.snp_meth.tsv
        col=()
    fi
done

for path in $(ls RCA/Normal/*.txt)
do
    tem=22

    # 提取名称
    m=${path/.RCA.transcript.txt/}
    name=$(basename ${m})

    col=()
    col+=("$tem")
    
    name2=$(awk '$1 == "'$name'" {print $2}' meth/Ala.nor.meth.txt)
    if [ -z "$name2" ]; then
        echo "$name meth is empty"
        for i in {1..2425}
        do
            col+=("")
        done
    else
        # 处理 CG, CHG, CHH 数据
        process_cx_coverage "CG" 262 "Normal_me/bismark/$name2/$name2.nameSorted.deduplicated.CG.CX_report.RCA.promoter_terminator.txt"
        process_cx_coverage "CHG" 337 "Normal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHG.CX_report.RCA.promoter_terminator.txt"
        process_cx_coverage "CHH" 1826 "Normal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHH.CX_report.RCA.promoter_terminator.txt"
    fi

    # 处理 snp 数据
    p=$(awk -v var="$name" -F',' '($1 == var || $11 == var) { print $30 }' Normal.csv | uniq)

    q=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk -F'[()]' '{print $2}')

    # snp_1 到 snp_104
    for snp in {1..104}
    do
        j=$(awk -v i="$snp" 'NR > 1 && NR == i + 1 {print $2}' snp/filtered_snp_matrix.csv)
        eval snp_$snp=$(awk -v q="$q" -v j="$j" 'NR==1 {for(i=1;i<=NF;i++) if($i==q) {col=i; break}; next} col && $2==j {print $col}' snp/filtered_snp_matrix.csv)
        if [ -z "$(eval echo '$snp_'$snp)" ]; then
            col+=("")
        else
            col+=($(eval echo '$snp_'$snp))
        fi
    done

    # 计算第四列第二行到第四行的值的总和
    # sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)
    sum=$(awk 'NR==5 {print $4}' $path)
    
    # 计算第四列第三行到第四行的值的总和
    sum2=$(awk 'NR>=3 && NR<=4 {sum += $4} END {print sum}' $path)

    # 计算第四列第二行的值
    value2=$(awk 'NR==2 {print $4}' $path)

    # 计算第四列第三行的值
    value3=$(awk 'NR==3 {print $4}' $path)

    # 计算第四列第四行的值
    value4=$(awk 'NR==4 {print $4}' $path)

    # 调试输出
    echo "Processing file: $path"
    echo "Values: sum=$sum, value2=$value2, value3=$value3, value4=$value4"

    # 确保 value2, value3, value4 和 sum 不是空的
    if [ -n "$value2" ] && [ -n "$value3" ] && [ -n "$value4" ] && [ -n "$sum" ]; then
        # 确保所有变量转换为普通浮点数格式
        sum=$(echo $sum | awk '{printf "%.10f\n", $1}')
        sum2=$(echo $sum2 | awk '{printf "%.10f\n", $1}')
        value2=$(echo $value2 | awk '{printf "%.10f\n", $1}')
        value3=$(echo $value3 | awk '{printf "%.10f\n", $1}')
        value4=$(echo $value4 | awk '{printf "%.10f\n", $1}')

        # 使用 bc 计算
        persum2=$(echo "scale=4; ($sum2/$sum)*100" | bc)
        per2=$(echo "scale=4; ($value2/$sum)*100" | bc)
        per3=$(echo "scale=4; ($value3/$sum)*100" | bc)
        per4=$(echo "scale=4; ($value4/$sum)*100" | bc)

        # 添加前导零
        persum2=$(printf "%.4f" $persum2)
        per2=$(printf "%.4f" $per2)
        per3=$(printf "%.4f" $per3)
        per4=$(printf "%.4f" $per4)
    else
        echo "Warning: Missing values in $path"
        per2=0
        per3=0
        per4=0
    fi
    col+=($sum $value2 $sum2 $value3 $value4 $per2 $persum2 $per3 $per4)

    echo -e "$(IFS=$'\t'; echo "${col[*]}")" >> RCA/Alt.weighted.snp_meth.tsv
    col=()
done
