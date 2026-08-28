#!/bin/bash
cd /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen
grep "^>" <(gunzip -c GCF_000001735.4_TAIR10.1_genomic.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat GCF_000001735.4_TAIR10.1_rna.fa.gz GCF_000001735.4_TAIR10.1_genomic.fa.gz > gentrome.fa.gz
salmon index -t gentrome.fa.gz -d decoys.txt -p 8 -i ../salmon_index
# salmon index -p 8 -t /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_rna.fa.gz -i /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/salmon_index

ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana
rna_fasta=$ref/refgen/GCF_000001735.4_TAIR10.1_rna.fa.gz
gff_file=$ref/refgen/GCF_000001735.4_TAIR10.1_genomic.gff

# 创建映射目录
mkdir -p $ref/mapping

# 提取rna.fa.gz文件中的转录本ID
echo "提取转录本ID..."
zcat $rna_fasta | grep "^>" | sed 's/^>//' | awk '{print $1}' > $ref/mapping/transcript_ids.txt

# 计算转录本ID数量
transcript_count=$(wc -l < $ref/mapping/transcript_ids.txt)
echo "从FASTA文件中提取到 $transcript_count 个转录本ID"

# 创建空的输出文件
> $ref/mapping/tx2gene.tsv

# 对每个转录本ID循环处理
echo "开始匹配转录本ID与基因ID..."
count=0
matched=0

while IFS= read -r transcript_id; do
    count=$((count+1))
    
    # 每处理100个ID打印一次进度
    if [ $((count % 100)) -eq 0 ]; then
        echo "已处理 $count / $transcript_count 个转录本ID，找到匹配: $matched"
    fi
    
    # 在GFF文件中查找该转录本ID
    # 使用grep匹配包含该转录本ID的行，然后用awk提取基因ID
    gene_id=$(grep -P "$transcript_id" $gff_file | grep -P "RNA\t" | 
              awk -F'\t' '{
                if($9 ~ /Parent=gene-/) {
                  match($9, /Parent=gene-([^;]+)/, arr);
                  print arr[1];
                  exit; # 只取第一个匹配结果
                }
              }')
    
    # 如果未找到基因ID，则尝试使用transcript模式
    if [ -z "$gene_id" ]; then
        gene_id=$(grep -P "$transcript_id" $gff_file | grep -P "transcript\t" | 
                  awk -F'\t' '{
                    if($9 ~ /Parent=gene-/) {
                      match($9, /Parent=gene-([^;]+)/, arr);
                      print arr[1];
                      exit; # 只取第一个匹配结果
                    }
                  }')
    fi
    
    # 如果找到了基因ID，则输出映射关系
    if [ ! -z "$gene_id" ]; then
        echo -e "$transcript_id\t$gene_id" >> $ref/mapping/tx2gene.tsv
        matched=$((matched+1))
    fi
done < $ref/mapping/transcript_ids.txt

# 输出最终映射文件的统计信息
echo "总共处理了 $count 个转录本ID"
echo "成功匹配了 $matched 个转录本ID到基因ID的映射关系"

# 检查生成的映射文件
echo "前10行映射文件内容:"
head $ref/mapping/tx2gene.tsv

# 统计信息
echo "映射文件总行数:"
wc -l $ref/mapping/tx2gene.tsv

# 清理临时文件
rm $ref/mapping/transcript_ids.txt

# 如果映射结果太少，添加警告
if [ $matched -lt 100 ]; then
    echo "警告: 找到的映射关系数量很少，可能需要检查数据格式或调整匹配方法"
    
    # 打印一些示例数据帮助调试
    echo "FASTA文件前3行示例:"
    zcat $rna_fasta | head -3
    
    echo "GFF文件前3行示例:"
    head -3 $gff_file
fi

# 无用代码示例
# ref=/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana
# # 创建映射目录
# mkdir -p $ref/mapping

# # 从GFF文件提取转录本和基因信息，针对指定的RefSeq格式
# zcat zcat /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_rna.fa.gz
# grep -P "RNA\t" $ref/refgen/GCF_000001735.4_TAIR10.1_genomic.gff | \
#   awk -F'\t' '{
#     info=$9;
#     transcript="";
#     gene_id="";
    
#     # 提取转录本ID - RefSeq格式: ID=rna-NM_001331245.1
#     if(match(info, /ID=rna-([^;]+)/, t)) {
#       transcript=t[1];
#     }

    
#     # 仅提取Parent=gene-格式的基因ID
#     if(match(info, /Parent=gene-([^;]+)/, p)) {
#       gene_id=p[1];
#     }
    
#     # 确保转录本ID和基因ID都不为空再输出
#     if(transcript != "" && gene_id != "") {
#       print transcript "\t" gene_id;
#     }
#   }' > $ref/mapping/tx2gene.tsv
# # 检查生成的映射文件
# echo "前10行映射文件内容:"
# head $ref/mapping/tx2gene.tsv
# echo "映射文件总行数:"
# wc -l $ref/mapping/tx2gene.tsv