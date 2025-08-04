#!/bin/bash

#for i in $(ls */bismark/*/*.CG.bw)
#do
    #file=$(basename $i)
    #file=$(echo $file | sed 's/\.nameSorted\.deduplicated//')
    #cp $i meth/$file
#done
#SRX1734567 16 Cvi-0 (6911) SRX248646
#SRX1734775 39.6589 IP-Men-2 (9556) SRX1664782
#SRX1734563 63.324 Bil-5 (6900) SRX1664577 

#cp Normal_me/bismark/SRX248646/*.CG.bw meth/low.bw
#cp Normal_me/bismark/SRX1664782/*.CG.bw meth/mid.bw
#cp Normal_me/bismark/SRX1664577/*.CG.bw meth/high.bw

#cp Normal_me/bismark/SRX248646/*.CHG.bw meth/low.bw
#cp Normal_me/bismark/SRX1664782/*.CHG.bw meth/mid.bw
#cp Normal_me/bismark/SRX1664577/*.CHG.bw meth/high.bw

cp Normal_me/bismark/SRX248646/*.CHH.bw meth/low.bw
cp Normal_me/bismark/SRX1664782/*.CHH.bw meth/mid.bw
cp Normal_me/bismark/SRX1664577/*.CHH.bw meth/high.bw

# 定义处理每个文件对或单文件的函数
#process_bw() {
    #local bigwig_files="$1"
    #f=$(basename $bigwig_files)
    #f=${f/.CG.bw/}

    computeMatrix scale-regions -S meth/low.bw meth/mid.bw meth/high.bw \
                                -R gene_RCA \
                                --numberOfProcessors 10 \
                                --beforeRegionStartLength 2000 \
                                --regionBodyLength 3000 \
                                --afterRegionStartLength 2000 \
                                --missingDataAsZero \
                                --binSize 100 \
                                -o meth/CHH.matrix.mat.gz
#}

# 导出函数以便 GNU Parallel 使用
#export -f process_bw


# 找出所有需要处理的文件

#parallel -j 128 process_bw ::: $(find meth -name *.CG.bw)

#total_files=$(find meth -name *.CG.matrix.mat.gz)
#合并所有的文件
#computeMatrixOperations rbind -m $total_files -o meth/CG.combined_matrix.mat.gz

# 绘制热图
plotHeatmap -m meth/CHH.matrix.mat.gz \
            --perGroup \
            --colorMap RdBu_r \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "RCA CHH_methylation level" \
            -out meth/CHH.png
			

#plotProfile -m meth/CG.matrix.mat.gz \
			#--perGroup \
			#--regionsLabel "RCA CG_methylation level" \
            #--outFileNameData meth/gene.profile.tab \
			#-out meth/CG.png

#plotProfile -m meth/CG.matrix.mat.gz \
            #--perGroup \
            #--regionsLabel "RCA CG_methylation level" \
            #--plotType heatmap \
            #-out meth/CG.headmap.png

#rm -rf meth/*CG.bw meth/*CG.matrix.mat.gz
rm -rf meth/low.bw meth/mid.bw meth/high.bw

