#!/bin/bash
for i in $(ls RCA3/*/*/*/*/*csv)
do
    m=$(dirname $i)
    if [ -d $m ]; then
        echo "dir exists"
    else
        mkdir $m
    fi
    for j in $(awk -F',' 'NR>1 {print $1}' $i; awk -F',' 'NR>1 {print $11}' $i)
    do
        #mv RCA/Abnormal/$j.RCA.transcript_count.txt $m
        mv RCA3/Abnormal/$j.RCA.transcript.txt $m
    done
done

for i in $(ls RCA3/*/*/*/*csv)
do
    m=$(dirname $i)
    if [ -d $m ]; then
        echo "dir exists"
    else
        mkdir $m
    fi
    for j in $(awk -F',' 'NR>1 {print $1}' $i; awk -F',' 'NR>1 {print $11}' $i)
    do
        mv RCA3/Normal/$j.RCA.transcript.txt $m
    done
done

for i in $(ls RCA3/*/*/*/*/*csv)
do
    rm -rf $i
done

for i in $(ls RCA3/*/*/*/*csv)
do
    rm -rf $i
done
exit 0