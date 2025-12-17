#!/bin/bash
# 初始化列名数组
columns=("tem")

# 动态生成 CG_1 到 CG_226
for i in {1..254}; do
    columns+=("CG_$i")
done

# 动态生成 CHG_1 到 CHG_301
for i in {1..325}; do
    columns+=("CHG_$i")
done

# 动态生成 CHH_1 到 CHH_1379
for i in {1..1774}; do
    columns+=("CHH_$i")
done

# 动态生成 promoter、terminator数据
# columns+=("CG_promoter" "CG" "CG_terminator" "CHG_promoter" "CHG" "CHG_terminator" "CHH_promoter" "CHH" "CHH_terminator")

columns+=("total" "α" "β" "β1" "β2" "α/%" "β/%" "β1/%" "β2/%")
# 写入表头
echo -e "$(IFS=$'\t'; echo "${columns[*]}")" > RCA/Alt.meth.tsv

# 定义函数处理 bedGraph 文件
process_bedgraph() {
    local prefix=$1  # CG, CHG, CHH
    local max=$2     # 最大迭代数
    local folder=$3  # 文件夹路径
    for i in $(seq 1 $max); do
        j=$(awk "NR == $i {print \$3}" meth/${prefix}.CV.filter.bedGraph)
        value=$(awk "\$3 == \"$j\" {print \$4}" "$folder")
        if [ -z "$value" ]; then
            value=""
        fi
        col+=("$value")
    done
}

# 定义函数处理 bedGraph1 文件
process_bedgraph1() {
    local prefix=$1  # CG, CHG, CHH
    local max=$2     # 最大迭代数
    local folder=$3  # 文件夹路径
    for i in $(seq 1 $max); do
        j=$(awk "NR == $i {print \$3}" meth/${prefix}.CV.filter.bedGraph)
        value=$(awk "\$3 == \"$j\" {print \$4}" "$folder")
        if [ -z "$value" ]; then
            value=""
        fi
        col1+=("$value")
    done
}

# 定义函数处理 bedGraph2 文件
process_bedgraph2() {
    local prefix=$1  # CG, CHG, CHH
    local max=$2     # 最大迭代数
    local folder=$3  # 文件夹路径
    for i in $(seq 1 $max); do
        j=$(awk "NR == $i {print \$3}" meth/${prefix}.CV.filter.bedGraph)
        value=$(awk "\$3 == \"$j\" {print \$4}" "$folder")
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
        process_bedgraph1 "CG" 254 "Abnormal_me/bismark/${name_array[0]}/${name_array[0]}.nameSorted.deduplicated.CG.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph1 "CHG" 325 "Abnormal_me/bismark/${name_array[0]}/${name_array[0]}.nameSorted.deduplicated.CHG.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph1 "CHH" 1774 "Abnormal_me/bismark/${name_array[0]}/${name_array[0]}.nameSorted.deduplicated.CHH.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph2 "CG" 254 "Abnormal_me/bismark/${name_array[1]}/${name_array[1]}.nameSorted.deduplicated.CG.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph2 "CHG" 325 "Abnormal_me/bismark/${name_array[1]}/${name_array[1]}.nameSorted.deduplicated.CHG.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph2 "CHH" 1774 "Abnormal_me/bismark/${name_array[1]}/${name_array[1]}.nameSorted.deduplicated.CHH.RCA.promoter_terminator.processed.bedGraph"

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

        echo -e "$(IFS=$'\t'; echo "${col1[*]}")" >> RCA/Alt.meth.tsv
        echo -e "$(IFS=$'\t'; echo "${col2[*]}")" >> RCA/Alt.meth.tsv
        col1=()
        col2=()

    else
        col+=("$tem")
        if [ -z "$name2" ]; then
            for i in {1..2353}
            do
                col+=("")
            done
        else
            # 处理 CG, CHG, CHH 数据
            process_bedgraph "CG" 254 "Abnormal_me/bismark/$name2/$name2.nameSorted.deduplicated.CG.RCA.promoter_terminator.processed.bedGraph"
            process_bedgraph "CHG" 325 "Abnormal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHG.RCA.promoter_terminator.processed.bedGraph"
            process_bedgraph "CHH" 1774 "Abnormal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHH.RCA.promoter_terminator.processed.bedGraph"
        fi

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

        echo -e "$(IFS=$'\t'; echo "${col[*]}")" >> RCA/Alt.meth.tsv
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
        for i in {1..2353}
        do
            col+=("")
        done
    else
        # 处理 CG, CHG, CHH 数据
        process_bedgraph "CG" 254 "Normal_me/bismark/$name2/$name2.nameSorted.deduplicated.CG.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph "CHG" 325 "Normal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHG.RCA.promoter_terminator.processed.bedGraph"
        process_bedgraph "CHH" 1774 "Normal_me/bismark/$name2/$name2.nameSorted.deduplicated.CHH.RCA.promoter_terminator.processed.bedGraph"
    fi

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

    echo -e "$(IFS=$'\t'; echo "${col[*]}")" >> RCA/Alt.meth.tsv
    col=()
done
