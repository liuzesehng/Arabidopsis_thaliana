#!/bin/bash
#RCA数据整理
echo -e "Tem/℃\tLat\tLong\talt/m\tTPM\tα_TPM/%\tβ_TPM/%\tβ1_TPM/%\tβ2_TPM/%" > Alt.RCA3.tsv
for path in $(ls Abnormal/*/*/*/*.txt)
do    
    # 使用cut命令和斜杠（/）作为分隔符，提取第2个字段
    tem=$(echo $path | cut -d'/' -f2)
    tem=${tem/C}

    # 提取名称
    m=${path/.RCA.transcript.txt/}
    name=$(basename ${m})

    # 使用awk选择第二列等于name的行的第四列
    latitude=$(awk '$1 == "'$name'" {print $2}' ../Ala.ab.location.txt)

    # 使用awk选择第二列等于name的行的第五列
    longitude=$(awk '$1 == "'$name'" {print $3}' ../Ala.ab.location.txt)

    # 使用awk选择第二列等于name的行的第五列
    alt=$(awk '$1 == "'$name'" {print $4}' ../Ala.ab.location.txt)

    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)

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
    echo -e "$tem\t$latitude\t$longitude\t$alt\t$sum\t$per2\t$persum2\t$per3\t$per4" >> Alt.RCA3.tsv
done

for path in $(ls Normal/*/*/*.txt)
do
    tem=22

    # 提取名称
    m=${path/.RCA.transcript.txt/}
    name=$(basename ${m})

    # 使用awk选择第二列等于name的行的第四列
    latitude=$(awk '$1 == "'$name'" {print $2}' ../Ala.normal.location.txt)

    # 使用awk选择第二列等于name的行的第五列
    longitude=$(awk '$1 == "'$name'" {print $3}' ../Ala.normal.location.txt)

    # 使用awk选择第二列等于name的行的第五列
    alt=$(awk '$1 == "'$name'" {print $4}' ../Ala.normal.location.txt)

    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)

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
    echo -e "$tem\t$latitude\t$longitude\t$alt\t$sum\t$per2\t$persum2\t$per3\t$per4" >> Alt.RCA3.tsv
done


exit 0

#RCA数据整理
echo -e "温度/℃\t纬度\t经度\tTPM\tα_TPM/%\tβ_TPM/%\tβ2_TPM/%" > RCA.tsv
for path in $(ls Abnormal/*/*/*/*.txt)
do
    # 使用cut命令和斜杠（/）作为分隔符，提取第4个字段
    tem=$(echo $path | cut -d'/' -f2)
    tem=${tem/C}
    country=$(echo $path | cut -d'/' -f4)

    # 使用awk选择第二列等于country的行的第四列
    latitude=$(awk -F, '$2 == "'$country'" {print $4}' ../country.csv | tr -d '\n' | tr -d '\r')

    # 使用awk选择第二列等于country的行的第五列
    longitude=$(awk -F, '$2 == "'$country'" {print $5}' ../country.csv | tr -d '\n' | tr -d '\r')

    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)

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
        per2=$(echo "scale=2; ($value2/$sum)*100" | bc)
        per3=$(echo "scale=2; ($value3/$sum)*100" | bc)
        per4=$(echo "scale=2; ($value4/$sum)*100" | bc)
    else
        echo "Warning: Missing values in $path"
        per2=0
        per3=0
        per4=0
    fi

    # 将所有数据合并成一行，并确保没有换行符
    echo -e "$tem\t$latitude\t$longitude\t$sum\t$per2\t$per3\t$per4" >> RCA.tsv
done

for path in $(ls Normal/*/*/*.txt)
do
    # 使用cut命令和斜杠（/）作为分隔符，提取第3个字段
    tem=22
    country=$(echo $path | cut -d'/' -f3)

    # 使用awk选择第二列等于country的行的第四列
    latitude=$(awk -F, '$2 == "'$country'" {print $4}' ../country.csv | tr -d '\n' | tr -d '\r')

    # 使用awk选择第二列等于country的行的第五列
    longitude=$(awk -F, '$2 == "'$country'" {print $5}' ../country.csv | tr -d '\n' | tr -d '\r')

    # 计算第四列第二行到第四行的值的总和
    sum=$(awk 'NR>=2 && NR<=4 {sum += $4} END {print sum}' $path)

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
        per2=$(echo "scale=2; ($value2/$sum)*100" | bc)
        per3=$(echo "scale=2; ($value3/$sum)*100" | bc)
        per4=$(echo "scale=2; ($value4/$sum)*100" | bc)
    else
        echo "Warning: Missing values in $path"
        per2=0
        per3=0
        per4=0
    fi


    # 将所有数据合并成一行，并确保没有换行符
    echo -e "$tem\t$latitude\t$longitude\t$sum\t$per2\t$per3\t$per4" >> RCA.tsv
done



