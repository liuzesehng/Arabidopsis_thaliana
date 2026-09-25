#！/bin/bash
for file in $(cat wheat.genome.txt | sed -n '26,34p'); do
    name=$(echo $file | cut -d'_' -f5)
    if [ -f $name.fasta.gz ] && [ ! -f $name.fasta.gz.aria2 ]; then
        echo "$name.fasta.gz exists, skip"
        continue
    fi
    # 使用 aria2 下载
    aria2c -x 16 -s 16 -k 1M \
        --file-allocation=none \
        --retry-wait=60 \
        --max-tries=10 \
        --timeout=600 \
        --connect-timeout=60 \
        --conditional-get=true \
        --allow-overwrite=true \
        --auto-file-renaming=false \
        --max-overall-download-limit=0 \
        --max-download-limit=0 \
        --continue=true \
        --console-log-level=notice \
        --summary-interval=1 \
        "$file" \
        -o "$name.fasta.gz"
done
