#!/bin/bash
#RCA数据整理
awk -F'\t' -v OFS='\t' 'NR==1{ $1="Temperature"; print $0, "total", "α", "β", "β1", "β2", "α/%", "β/%", "β1/%", "β2/%"; exit }' Ala.ab.location.climate.txt > RCA/RCA.climate.tsv
for path in $(ls RCA/Abnormal/*/*.txt)
do    
    # 使用cut命令和斜杠（/）作为分隔符，提取第2个字段
    tem=$(echo $path | cut -d'/' -f3)
    tem=${tem/C}

    # 提取名称
    m=${path/.RCA.transcript.txt/}
    name=$(basename ${m})

    # 获取匹配样本的位置信息（第2列及以后）
    location_fields=$(awk -F'\t' -v OFS='\t' -v name="$name" '$1==name { $1=""; sub(/^\t/, ""); print; exit }' Ala.ab.location.climate.txt)
    if [ -z "$location_fields" ]; then
        echo "Warning: Missing location fields for $name in Ala.ab.location.climate.txt"
        location_fields="\t\t"  # 用适当数量的制表符填充以保持列对齐
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

    # 将所有数据合并成一行，并确保没有换行符
    echo -e "$tem\t$location_fields\t$sum\t$value2\t$sum2\t$value3\t$value4\t$per2\t$persum2\t$per3\t$per4" >> RCA/RCA.climate.tsv
done

for path in $(ls RCA/Normal/*.txt)
do
    tem=22

    # 提取名称
    m=${path/.RCA.transcript.txt/}
    name=$(basename ${m})

    # 获取匹配样本的位置信息（第2列及以后）
    location_fields=$(awk -F'\t' -v OFS='\t' -v name="$name" '$1==name { $1=""; sub(/^\t/, ""); print; exit }' Ala.normal.location.climate.txt)
    if [ -z "$location_fields" ]; then
        echo "Warning: Missing location fields for $name in Ala.normal.location.climate.txt"
        location_fields="\t\t"  # 用适当数量的制表符填充以保持列对齐
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


    # 将所有数据合并成一行，并确保没有换行符
    echo -e "$tem\t$location_fields\t$sum\t$value2\t$sum2\t$value3\t$value4\t$per2\t$persum2\t$per3\t$per4" >> RCA/RCA.climate.tsv
done
