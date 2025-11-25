#!/bin/bash
for i in $(ls salmon/Normal)
do
    p=$(awk -v var="$i" -F, '($1 == var || $11 == var) { print $30 }' Normal.csv | uniq)

    q=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk -F'[()]' '{print $2}')

    o=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk '{print $1}' | sed 's/"//')

    output=$(awk -v var="$q" -v new_var="$i" -F'\t' '($1 == var) { print new_var "\t" $3 "\t" $4 "\t" $5 }' species.location.modified.processed.txt)

    if [ -z "$output" ]; then
        echo "$i" >> Ala.normal.location.txt
        echo "$q $i $o not exist"
    else
        echo "$output" >> Ala.normal.location.txt
    fi
done
sed -i '1iSample\tlat\tlong\talt' Ala.normal.location.txt

# exit 0

for i in $(ls salmon/Abnormal)
do
    p=$(awk -v var="$i" -F, '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)

    q=$(cat biosample_result.txt | grep -A 4 $p | sed -n '5,5p' | sed 's/.*\/accession number="\([^"]*\)".*/\1/')

    o=$(cat biosample_result.txt | grep -A 5 $p | sed -n '6,6p' | sed 's/.*\/accession name="\([^"]*\)".*/\1/')

    output=$(awk -v var="$q" -v new_var="$i" -F'\t' '($1 == var) { print new_var "\t" $3 "\t" $4 "\t" $5 }' species.location.modified.processed.txt)
    
    if [ -z "$output" ]; then
        echo "$i" >> Ala.ab.location.txt
        echo "$q $i $o not exist"
    else
        echo "$output" >> Ala.ab.location.txt
    fi
done
# exit 0
sed -i '1iSample\tlat\tlong\talt' Ala.ab.location.txt






