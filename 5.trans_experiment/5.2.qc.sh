#!/bin/bash
process_fastp(){
    local file="$1"
    dir=$(dirname "$file")
    name=$(basename "$file" .fq.gz)
    file2=${file/_1.fq.gz/_2.fq.gz}
    name2=$(basename "$file2" .fq.gz)
    if [ ! -d $dir ]; then
        mkdir -p $dir
    fi

    fit=$dir/$name".fit.fastq.gz"
    fit2=$dir/$name2".fit.fastq.gz"

    fastp -z 9 -w 4 -i $file -I $file2 -o $fit -O $fit2 \
        --disable_quality_filtering \
        --disable_length_filtering \
        --html $dir/$name"_fastp.html" \
        --json $dir/$name"_fastp.json"


}
export -f process_fastp
ls rnadata/*/*fq.gz | grep -v "_2.fq.gz" | parallel -j 16 process_fastp {}

