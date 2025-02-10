#!/bin/bash
p+=(0)
p+=(1)
for i in $(ls Normal_me/*sra)
do
    n=${i/.sra/}
    n1=${i/\/*/}
    n2=${n/*\//}
    n3=$(awk -v var="$n2" -F',' '$1 == var {print $16}' Normal_me/Normal.me.csv)
    n4=${m/*\//}
    n5=$(awk -v var="$n4" -F',' '$1 == var {print $16}' Normal_me/Normal.me.csv)
    if [ "$n3" = "PAIRED" ]; then
        fastq-dump --split-files $i -O $n1

        j=$(awk -v var="$n2" -F',' '$1 == var {print $11}' Normal_me/Normal.me.csv)
        if [ "$j" = "${p[-1]}" ]; then
            cat "$m"_1.fastq >> $n1/"$j"_1.fastq
            cat "$m"_2.fastq >> $n1/"$j"_2.fastq
            rm -rf "$m"_1.fastq "$m"_2.fastq
        else
            if [ "${p[-1]}" = "${p[-2]}" ]; then
                if [ "$n5" = "PAIRED" ]; then
                    cat "$m"_1.fastq >> $n1/${p[-1]}_1.fastq
                    cat "$m"_2.fastq >> $n1/${p[-1]}_2.fastq
                    rm -rf "$m"_1.fastq "$m"_2.fastq
                else
                    cat "$m".fastq >> $n1/${p[-1]}.fastq
                    rm -rf "$m".fastq
                fi    
            fi
        fi
    else
        fastq-dump $i -O $n1

        j=$(awk -v var="$n2" -F',' '$1 == var {print $11}' Normal_me/Normal.me.csv)

        if [ "$j" = "${p[-1]}" ]; then
            cat $m.fastq >> $n1/$j.fastq
            rm -rf $m.fastq
        else
            if [ "${p[-1]}" = "${p[-2]}" ]; then
                if [ "$n5" = "PAIRED" ]; then
                    cat "$m"_1.fastq >> $n1/${p[-1]}_1.fastq
                    cat "$m"_2.fastq >> $n1/${p[-1]}_2.fastq
                    rm -rf "$m"_1.fastq "$m"_2.fastq
                else
                    cat "$m".fastq >> $n1/${p[-1]}.fastq
                    rm -rf "$m".fastq
                fi
            fi
        fi
    fi
    m=$n
    p+=($j) 
done
cat "$m".fastq >> $n1/${p[-1]}.fastq
#cat "$m"_2.fastq >> $n1/${p[-1]}_2.fastq
rm -rf "$m".fastq

exit 0

# Download the data from the SRA database
for i in $(ls  Normal_me/*csv)
do
    n=${i%/*}
    if [ -d ./$n ]; then
        echo "dir exists"
    else
        mkdir ./$n
    fi
    for j in $(awk -F',' 'NR>1 {print $10}' $i)
    do
        m=${j##*/}
        p=${m%%.*}
        #q=$(awk -v var="$p" -F',' '$1 == var {print $11}' $i)
        if [ -f $n/$p.sra ]; then
            echo "file exists"
        else
            # 尝试下载，直到成功为止
            while true; do
                wget -c "$j"
                if [ $? -eq 0 ]; then
                    echo "下载成功"
                    break
                else
                    echo "下载失败，重试中..."
                    sleep 5  # 等待 5 秒后重试
                fi
            done
            mv $m $n/$p.sra
        fi
    done
done

for i in $(ls Normal_me/*csv)
do
    n=${i%/*}
    if [ -d ./$n ]; then
        echo "dir exists"
    else
        mkdir ./$n
    fi
    for j in $(awk -F',' 'NR>1 {print $10}' $i)
    do
        m=${j##*/}
        p=${m%%.*}
        #q=$(awk -v var="$p" -F',' '$1 == var {print $11}' $i)
        if [ -f $n/$p.sra ]; then
            echo "file exists"
        else
            # 尝试下载，直到成功为止
            while true; do
                wget -c "$j"
                if [ $? -eq 0 ]; then
                    echo "下载成功"
                    break
                else
                    echo "下载失败，重试中..."
                    sleep 5  # 等待 5 秒后重试
                fi
            done
            mv $m $n/$p.sra
        fi
    done
done