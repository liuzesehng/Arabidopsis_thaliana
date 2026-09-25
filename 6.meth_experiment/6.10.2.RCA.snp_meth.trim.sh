#!/bin/bash
# 初始化列名数组
columns=("Temperature")

# 表头和取值函数使用同一参考文件的有效记录顺序，不写死位点数量。
for prefix in CG CHG CHH; do
    reference="meth_experiment/${prefix}.CV.before.numbered.tsv"
    if [[ ! -r "$reference" ]]; then
        echo "Error: Cannot read $reference" >&2
        exit 1
    fi
    # 第一列为已有位点编号，保留非连续编号及原始顺序。
    identifiers=$(awk 'NF && $1 !~ /^#/ {print $1}' "$reference") || exit 1
    while IFS= read -r identifier; do
        [[ -n "$identifier" ]] && columns+=("$identifier")
    done <<< "$identifiers"
done

columns+=("total" "α" "β" "α/%" "β/%")
# 按参考文件顺序匹配染色体和位点；缺失位点保留空单元格。
process_bedgraph() {
    local prefix=$1
    local folder=$2
    local values value
    if [[ ! -r "$folder" ]]; then
        echo "Error: Cannot read $folder" >&2
        return 1
    fi
    values=$(awk '
        NF && $1 !~ /^#/ && $1 != "track" && $1 != "browser" {
            if (FILENAME == ARGV[1]) {
                key = $1 SUBSEP $3
                if (!(key in data)) data[key] = $4
            } else {
                # numbered.tsv：第 1 列编号，第 2 列染色体，第 4 列终点。
                key = $2 SUBSEP $4
                # 标记前缀避免命令替换丢失末尾空值。
                print "value:" ((key in data) ? data[key] : "")
            }
        }
    ' "$folder" "meth_experiment/${prefix}.CV.before.numbered.tsv") || return 1
    while IFS= read -r value; do
        [[ -n "$value" ]] && col+=("${value#value:}")
    done <<< "$values"
    return 0
}

# 同温度内分别按样本编号数值排序，再按位置配对。
# 以 RNA 为输出基准，多出的 RNA 行保留表达值，甲基化列留空。
collect_samples() {
    local kind=$1 path sample temperature replicate
    local -a files=()
    if [[ "$kind" == meth ]]; then
        files=(methdata/*/*.R1.fit.fastq.gz)
    else
        files=(RCA/rnadata/*.RCA.transcript.txt RCA/rnadata/*.RCA.transition.txt)
    fi
    for path in "${files[@]}"; do
        sample=${path##*/}
        sample=${sample%%.*}
        sample=${sample//_/-}
        if [[ ! "$sample" =~ ^T([0-9]+)-([0-9]+)$ ]]; then
            echo "Error: Invalid sample name: $path" >&2
            return 1
        fi
        temperature=$((10#${BASH_REMATCH[1]}))
        replicate=$((10#${BASH_REMATCH[2]}))
        printf '%s\t%s\t%s\t%s\n' "$temperature" "$replicate" "$sample" "$path"
    done
}

shopt -s nullglob
meth_records=$(collect_samples meth) || exit 1
rna_records=$(collect_samples rna) || exit 1
if [[ -z "$rna_records" ]]; then
    echo "Error: No RNA input files in RCA/rnadata" >&2
    exit 1
fi
meth_records=$(printf '%s\n' "$meth_records" | sort -t $'\t' -k1,1n -k2,2n)
rna_records=$(printf '%s\n' "$rna_records" | sort -t $'\t' -k1,1n -k2,2n)
declare -A meth_counts=() rna_counts=() meth_names=() seen=()
while IFS=$'\t' read -r temperature replicate sample path; do
    [[ -n "$temperature" ]] || continue
    key="$temperature:$replicate"
    if [[ -n "${seen[$key]:-}" ]]; then
        echo "Error: Duplicate methylation sample: $sample" >&2
        exit 1
    fi
    seen[$key]=1
    rank=$((${meth_counts[$temperature]:-0} + 1))
    meth_counts[$temperature]=$rank
    meth_names["$temperature:$rank"]=$sample
done <<< "$meth_records"
seen=()
while IFS=$'\t' read -r temperature replicate sample path; do
    key="$temperature:$replicate"
    if [[ -n "${seen[$key]:-}" ]]; then
        echo "Error: Duplicate RNA sample: $sample" >&2
        exit 1
    fi
    seen[$key]=1
    rna_counts[$temperature]=$((${rna_counts[$temperature]:-0} + 1))
done <<< "$rna_records"
for temperature in "${!meth_counts[@]}"; do
    if ((${meth_counts[$temperature]} > ${rna_counts[$temperature]:-0})); then
        echo "Warning: T$temperature has extra methylation samples; only RNA-matched ranks will be output" >&2
    fi
done

# 写入表头
echo -e "$(IFS=$'\t'; echo "${columns[*]}")" > meth_experiment/Alt.all.snp_meth.tsv

# RNA 样本顺序即输出行顺序；列名 Temperature 保持兼容。
rna_counts=()
while IFS=$'\t' read -r temperature replicate rna_sample rna_file; do
    rank=$((${rna_counts[$temperature]:-0} + 1))
    rna_counts[$temperature]=$rank
    name=${meth_names["$temperature:$rank"]:-}
    col=("$rna_sample")
    if [[ -n "$name" ]]; then
        echo "RNA mapping: $name -> $rna_sample ($rna_file)"
        process_bedgraph "CG" "methdata/bismark/$name/$name.nameSorted.deduplicated.CG.RCA.promoter_terminator.processed.bedGraph" || exit 1
        process_bedgraph "CHG" "methdata/bismark/$name/$name.nameSorted.deduplicated.CHG.RCA.promoter_terminator.processed.bedGraph" || exit 1
        process_bedgraph "CHH" "methdata/bismark/$name/$name.nameSorted.deduplicated.CHH.RCA.promoter_terminator.processed.bedGraph" || exit 1
    else
        echo "RNA mapping: no methylation sample -> $rna_sample ($rna_file)"
        for ((i=1; i<=${#columns[@]}-6; i++)); do
            col+=("")
        done
    fi

    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR==5 {print $4}' "$rna_file")

    # 计算第四列第二行的值
    value2=$(awk 'NR==2 {print $4}' "$rna_file")

    # 计算第四列第三行的值
    value3=$(awk 'NR==3 {print $4}' "$rna_file")

    # 计算第四列第四行的值
    value4=$(awk 'NR==4 {print $4}' "$rna_file")

    # 计算第四列第三行到第四行的值的总和
    sum2=$(awk 'NR>=3 && NR<=4 {sum += $4} END {print sum}' "$rna_file")

    # 调试输出
    echo "Processing file: $rna_file"
    echo "Values: sum=$sum, value2=$value2, value3=$value3, value4=$value4"

    # 确保 value2, value3, value4 和 sum 不是空的
    if [ -n "$value2" ] && [ -n "$value3" ] && [ -n "$value4" ] && [ -n "$sum" ]; then
        # 确保所有变量转换为普通浮点数格式
        sum=$(echo $sum | awk '{printf "%.10f\n", $1}')
        sum2=$(echo $sum2 | awk '{printf "%.10f\n", $1}')
        value2=$(echo $value2 | awk '{printf "%.10f\n", $1}')
        value3=$(echo $value3 | awk '{printf "%.10f\n", $1}')
        value4=$(echo $value4 | awk '{printf "%.10f\n", $1}')

        if ! awk -v total="$sum" 'BEGIN {exit !(total > 0)}'; then
            echo "Error: Non-positive total in $rna_file" >&2
            exit 1
        fi

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
        echo "Error: Missing expression values in $rna_file" >&2
        exit 1
    fi
    col+=("$sum" "$value2" "$sum2" "$per2" "$persum2")

    echo -e "$(IFS=$'\t'; echo "${col[*]}")" >> meth_experiment/Alt.all.snp_meth.tsv
    col=()
done <<< "$rna_records"
