#!/bin/bash
process_fastp(){
    local file="$1"
    dir=$(dirname "$file")
    name=$(basename "$file" .fastq.gz)
    if [ "$name" = "SRR3465232_1" ]; then
        name=${name/_1/}
    fi
    file2=${file/_1/_2}
    if [ ! -d $dir/$name ]; then
        mkdir -p $dir/$name
    fi
    fit=$dir/$name/$name".fit.fastq.gz"
    fit1=$dir/$name/$name"_1.fit.fastq.gz"
    fit2=$dir/$name/$name"_2.fit.fastq.gz"

    if [ "$name" = "SRR3465232" ]; then
        fastp -z 9 -w 4 -i $file -I $file2 -o $fit1 -O $fit2 \
             --disable_quality_filtering \
             --disable_length_filtering \
             --html $dir/$name/$name"_fastp.html" \
             --json $dir/$name/$name"_fastp.json"
        # fastp -z 9 -w 4 -i $file -I $file2 -o $fit1 -O $fit2 -j $dir/$name/$name"_fastp.json" -h $dir/$name/$name"_fastp.html"
        # rm -rf $file $file2
    else
        fastp -z 9 -w 4 -i $file -o $fit \
             --disable_quality_filtering \
             --disable_length_filtering \
             --html $dir/$name/$name"_fastp.html" \
             --json $dir/$name/$name"_fastp.json"
        # fastp -z 9 -w 4 -i $file -o $fit -j $dir/$name/$name"_fastp.json" -h $dir/$name/$name"_fastp.html"
        # rm -rf $file
    fi
}
export -f process_fastp
ls Abnormal/*fastq.gz | grep -v "_2" | parallel -j 16 process_fastp {}
ls Normal/*fastq.gz | grep -v "_2" | parallel -j 16 process_fastp {}
# for i in $(ls Normal/SRR3465176.fastq.gz)
# do          
#     m=${i/*\//}
#     n=${m/.*gz/}
#     #if [ -d ./$n ]; then
#         #$echo "dir exists"
#     #else
#         #mkdir ./$n
#     #fi
#     fastqc -t 8 -o fastqc/Normal $i &
# done
wait
