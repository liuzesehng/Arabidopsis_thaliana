#!/bin/bash

# # 1.Download the data from the SRA database
# for i in $(ls Abnormal/*csv)
# do
#     n=${i%/*}
#     if [ -d ./$n ]; then
#         echo "dir exists"
#     else
#         mkdir ./$n
#     fi
#     for j in $(awk -F',' 'NR>1 {print $10}' $i)
#     do
#         m=${j##*/}
#         p=${m%%.*}
#         #q=$(awk -v var="$p" -F',' '$1 == var {print $11}' $i)
#         if [ -f $n/$p.sra ]; then
#             echo "file exists"
#         else
#             wget -c $j
#             mv $m $n/$p.sra
#         fi
#     done
# done
# exit 0


# for i in $(ls Normal/*csv)
# do
#     n=${i%/*}
#     if [ -d ./$n ]; then
#         echo "dir exists"
#     else
#         mkdir ./$n
#     fi
#     for j in $(awk -F',' 'NR>1 {print $10}' $i)
#     do
#         m=${j##*/}
#         p=${m%%.*}
#         #q=$(awk -v var="$p" -F',' '$1 == var {print $11}' $i)
#         if [ -f $n/$p.sra ]; then
#             echo "file exists"
#         else
#             wget -c $j
#             mv $m $n/$p.sra
#         fi
#     done
# done
# exit 0


# # 2.combind the data
# p+=(0)
# p+=(1)
# for i in $(ls Abnormal/*sra)
# do
#     n=${i/.sra/}
#     n1=${i/\/*/}
#     n2=${n/*\//}
#     if [ "$n2" = "SRR3465232" ]; then
#         fastq-dump --split-files $i -O $n1
#     else
#         fastq-dump $i -O $n1
#     fi
#     j=$(awk -v var="$n2" -F',' '$1 == var {print $11}' Abnormal/Abnormal.csv)
#     if [ "$j" = "${p[-1]}" ]; then
#         cat $m.fastq >> $n1/$j.fastq
#         rm -rf $m.fastq
#     else
#         if [ "${p[-1]}" = "${p[-2]}" ]; then    
#             cat $m.fastq >> $n1/${p[-1]}.fastq
#             rm -rf $m.fastq
#         fi
#     fi
#     m=$n
#     p+=($j) 
# done
# #rm -rf $m.fastq
# exit 0

# p+=(0)
# p+=(1)
# for i in $(ls Normal/*sra)
# do
#     n=${i/.sra/}
#     n1=${i/\/*/}
#     n2=${n/*\//}
#     if [ "$n2" = "SRR3465232" ]; then
#         fastq-dump --split-files $i -O $n1
#     else
#         fastq-dump $i -O $n1
#     fi
#     j=$(awk -v var="$n2" -F',' '$1 == var {print $11}' Normal/Normal.csv)
#     if [ "$j" = "${p[-1]}" ]; then
#         cat $m.fastq >> $n1/$j.fastq
#         rm -rf $m.fastq
#     else
#         if [ "${p[-1]}" = "${p[-2]}" ]; then    
#             cat $m.fastq >> $n1/${p[-1]}.fastq
#             rm -rf $m.fastq
#         fi
#     fi
#     m=$n
#     p+=($j) 
# done
# #rm -rf $m.fastq
# exit 0


# 3.compress the data
process_gzip(){
    local file=$1
    if [ -f "$file" ] && [ ! -f "$file.gz" ]; then
        echo "file exists: $file"
        gzip -9 "$file"
    else
        echo "file not exists or already compressed, skip: $file"
        return
    fi
}
export -f process_gzip

# ls Abnormal/*fastq | parallel -j 64 process_gzip {}
ls Normal/*fastq | parallel -j 64 process_gzip {}

exit 0




