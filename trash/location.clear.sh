#!/bin/bash


# 提取第4、5列和在剩余的列中提取包含数字的列
#awk '{fourth_column=$4; fifth_column=$5; numeric_columns=""; for (i=6; i<=NF; i++) {if ($i ~ /[0-9]/) {numeric_columns=numeric_columns $i " "}} print fourth_column, fifth_column, numeric_columns}' out.species.txt > processed_output.txt

# 使用 awk 删除只有空格的行，并将空格替换为制表符
#sed 's/[[:space:]]*$//' processed_output.txt | sed 's/ \+/\t/g' > species.location.txt

#rm -rf processed_output.txt

# 使用 awk 删除第3列和第4列的N、S、E、W，并将其转换为数字
#awk -F'\t' '{ if ($3 ~ /N|S/) { gsub(/N|S/, "", $3); split($3, arr, "-"); if (length(arr) == 2) { $3 = (arr[1] + arr[2]) / 2 } } if ($4 ~ /E|W/) { gsub(/E|W/, "", $4); split($4, arr, "-"); if (length(arr) == 2) { $4 = (arr[1] + arr[2]) / 2 } } print $0 }' OFS='\t' species.location.modified.txt > species.location.modified.processed.txt

#grep -A 1 "name" kml.txt | grep -B 1 "description" | grep -v '^--$' > kml.processed.txt 

#awk -F'[<>]' '/<name>/ { name = $3 } /<description>/ { n = split($0, arr, "><"); for (i = n; i > 0; i--) { if (arr[i] ~ /[0-9]/) { description = arr[i]; break } } print name, description }' kml.processed.txt > kml.final.txt

