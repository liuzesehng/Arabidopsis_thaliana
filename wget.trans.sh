#!/bin/bash
for i in $(ls Normal/*fastq)
do
    gzip -9 $i
    rm -rf $i
done
exit 0

p+=(0)
p+=(1)
for i in $(ls Normal/*sra)
do
    n=${i/.sra/}
    n1=${i/\/*/}
    n2=${n/*\//}
    if [ "$n2" = "SRR3465232" ]; then
        fastq-dump --split-files $i -O $n1
    else
        fastq-dump $i -O $n1
    fi
    j=$(awk -v var="$n2" -F',' '$1 == var {print $11}' Normal/Normal.csv)
    if [ "$j" = "${p[-1]}" ]; then
        cat $m.fastq >> $n1/$j.fastq
        rm -rf $m.fastq
    else
        if [ "${p[-1]}" = "${p[-2]}" ]; then    
            cat $m.fastq >> $n1/${p[-1]}.fastq
            rm -rf $m.fastq
        fi
    fi
    m=$n
    p+=($j) 
done
#rm -rf $m.fastq
exit 0

# Download the data from the SRA database
for i in $(ls Normal/*csv)
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
            wget -c $j
            mv $m $n/$p.sra
        fi
    done
done
exit 0



for i in $(ls Normal/*lite.sra)
do
    m=${i/.sralite/}
    n=${m/*\//}
    mv $i Normal/$n
done
exit 0

