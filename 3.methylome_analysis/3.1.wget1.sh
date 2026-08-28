#!/bin/bash

# 使用 GNU Parallel 并行压缩
ls Normal_me/*fastq | parallel -j 64 gzip -9

exit 0

# 分割文件
split -n l/10 fastq.txt segment

# 对每一份中的每个路径的文件进行压缩
for file in segment*; do
    while IFS= read -r line; do
        gzip -9 "$line" &
    done < "$file"
done
wait
rm -rf segment*
