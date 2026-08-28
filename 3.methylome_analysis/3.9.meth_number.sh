#!/bin/bash
for i in $(ls salmon/Abnormal)
do
    p=$(awk -v var="$i" -F',' '($1 == var || $11 == var) { print $30 }' Unnormal.csv | uniq)

    q=$(cat biosample_result.txt | grep -A 9 $p | grep "accession number=" | sed 's/.*\/accession number="\([^"]*\)".*/\1/')

    o=$(cat biosample_result.txt | grep -A 9 $p | grep "growth temperature=" | sed 's/.*\/growth temperature="\([^"]*\)".*/\1/')

    r=$(cat meth/biosample_result_me.txt | grep -B 5 "/accession number=\"$q\"" | grep "_rep")

    if [ -z "$r" ]; then
        m=$(cat meth/biosample_result_me.txt | grep -A 5 "/accession number=\"$q\"" | grep -A 3 "/growth temperature=\"$o\"" | grep "Accession:" | sed 's/.*Accession: \([^ ID]*\).*/\1/')

        output1=$(awk -v var=$m -F',' '($26 == var) { print $1 "\t" $11 }' Abnormal_me/Abnormal.me.csv | awk 'NR == 1 { first_col = $1 } NR == 2 { print $2; exit } END { if (NR < 2) print first_col }')

        if [ -z "$output1" ]; then
            echo -e "$p $q $o not exist"
            echo -e "$i" >> meth/Ala.abnor.meth.txt
        else
            echo -e "$i\t$output1\t" >> meth/Ala.abnor.meth.txt
        fi

    else
        m=$(cat meth/biosample_result_me.txt | grep -A 5 "/accession number=\"$q\"" | grep -A 3 "/growth temperature=\"$o\"" | grep "Accession:" | sed 's/.*Accession: \([^ ID]*\).*/\1/')
        m1=$(awk 'NR==1' <<< $m)
        m2=$(awk 'NR==2' <<< $m)

        output2=$(awk -v var=$m1 -F',' '($26 == var) { print $1 "\t" $11 }' Abnormal_me/Abnormal.me.csv | awk 'NR == 1 { first_col = $1 } NR == 2 { print $2; exit } END { if (NR < 2) print first_col }')
        output3=$(awk -v var=$m2 -F',' '($26 == var) { print $1 "\t" $11 }' Abnormal_me/Abnormal.me.csv | awk 'NR == 1 { first_col = $1 } NR == 2 { print $2; exit } END { if (NR < 2) print first_col }')

        if [ -z "$output2" ]; then
            echo -e "$p $q $o not exist"
            echo -e "$i" >> meth/Ala.abnor.meth.txt
        else
            echo -e "$i\t$output2\t" >> meth/Ala.abnor.meth.txt
        fi

        if [ -z "$output3" ]; then
            echo -e "$p $q $o not exist"
            echo -e "$i" >> meth/Ala.abnor.meth.txt
        else
            echo -e "$i\t$output3\t" >> meth/Ala.abnor.meth.txt
        fi

    fi

done

for i in $(ls salmon/Normal)
do

    p=$(awk -v var="$i" -F',' '($1 == var || $11 == var) { print $30 }' Normal.csv | uniq)

    q=$(awk -v var="$p" -F',' '($1 == var) {print $2}' sample.csv | awk -F'[()]' '{print $2}')

    m=$(cat meth/biosample_result_nor_me.txt | grep -A 3 "/ecotype id=\"$q\"" | grep "Accession:" | sed 's/.*Accession: \([^ ID]*\).*/\1/')

    output=$(awk -v var=$m -F',' '($26 == var) { print $1 "\t" $11 }' Normal_me/Normal.me.csv | awk 'NR == 1 { first_col = $1 } NR == 2 { print $2; exit } END { if (NR < 2) print first_col }')

    if [ -z "$output" ]; then
        echo -e "$p $q not exist"
        echo -e "$i" >> meth/Ala.nor.meth.txt
    else
        echo -e "$i\t$output\t" >> meth/Ala.nor.meth.txt
    fi
done
