for i in $(ls Normal/SRR3465176.fastq.gz)
do          
    m=${i/*\//}
    n=${m/.*gz/}
    #if [ -d ./$n ]; then
        #$echo "dir exists"
    #else
        #mkdir ./$n
    #fi
    fastqc -t 8 -o fastqc/Normal $i &
done
wait