#!/bin/bash
source ~/.bashrc
conda activate bismark
for i in CG CHG CHH
do
    cut -d $'\t' -f 1-3,5 meth/$i.CV.bedGraph > meth/$i.bedGraph
    bedGraphToBigWig meth/$i.bedGraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt meth/$i.bw
    rm -rf meth/$i.bedGraph
done
conda activate deeptools
computeMatrix scale-regions -S meth/CG.bw meth/CHG.bw meth/CHH.bw \
                            -R gene_RCA \
                            --numberOfProcessors 10 \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --binSize 100 \
                            -o meth/meth.matrix.mat.gz

plotHeatmap -m meth/meth.matrix.mat.gz \
            --perGroup \
            --colorMap RdBu_r \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "Rca methylation"  \
            -out meth/meth.cv.pdf
rm -rf meth/CG.bw meth/CHG.bw meth/CHH.bw meth/meth.matrix.mat.gz

conda activate bismark
for i in CG CHG CHH
do
    cut -d $'\t' -f 1-3,4 meth/$i.CV.bedGraph > meth/$i.bedGraph
    bedGraphToBigWig meth/$i.bedGraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt meth/$i.bw
    rm -rf meth/$i.bedGraph
done

conda activate deeptools

computeMatrix scale-regions -S meth/CG.bw meth/CHG.bw meth/CHH.bw \
                            -R gene_RCA \
                            --numberOfProcessors 10 \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --binSize 100 \
                            -o meth/meth.matrix.mat.gz

plotHeatmap -m meth/meth.matrix.mat.gz \
            --perGroup \
            --colorMap RdBu_r \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "Rca methylation" \
            -out meth/meth.pdf
rm -rf meth/CG.bw meth/CHG.bw meth/CHH.bw

exit 0

computeMatrix scale-regions -S meth/CG.bw meth/CHG.bw meth/CHH.bw \
                            -R gene_RCA \
                            --numberOfProcessors 10 \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --binSize 100 \
                            -o meth/meth.matrix.mat.gz

# 绘制热图
plotHeatmap -m meth/meth.matrix.mat.gz \
            --perGroup \
            --whatToShow 'heatmap and colorbar' \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "RCA methylation level" \
            --dpi 1200 \
            -out meth/meth.heat.mean.pdf
			

plotProfile -m meth/meth.matrix.mat.gz \
			--perGroup \
            --plotHeight 7 \
            --plotWidth 15 \
			--regionsLabel "RCA methylation level" \
            --dpi 1200 \
			-out meth/meth.mean.pdf
#rm -rf meth/*CG.bw meth/*CG.matrix.mat.gz
rm -rf meth/CG.bw meth/CHG.bw meth/CHH.bw
exit 0

for i in CG CHG CHH
do
    cut -d $'\t' -f 1-3,4 meth/$i.CV.bedGraph > meth/$i.bedGraph
    bedGraphToBigWig meth/$i.bedGraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt meth/$i.bw
    rm -rf meth/$i.bedGraph
done
exit 0

computeMatrix scale-regions -S meth/CG.bw meth/CHG.bw meth/CHH.bw \
                            -R gene_RCA \
                            --numberOfProcessors 10 \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --binSize 100 \
                            -o meth/meth.matrix.mat.gz

# 绘制热图
plotHeatmap -m meth/meth.matrix.mat.gz \
            --perGroup \
            --whatToShow 'heatmap and colorbar' \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "RCA methylation level" \
            --dpi 1200 \
            -out meth/meth.heat.pdf
			

plotProfile -m meth/meth.matrix.mat.gz \
			--perGroup \
            --plotHeight 7 \
            --plotWidth 15 \
			--regionsLabel "RCA methylation level" \
            --dpi 1200 \
			-out meth/meth.pdf

#plotProfile -m meth/CG.matrix.mat.gz \
            #--perGroup \
            #--regionsLabel "RCA CG_methylation level" \
            #--plotType heatmap \
            #-out meth/CG.headmap.png

#rm -rf meth/*CG.bw meth/*CG.matrix.mat.gz
rm -rf meth/CG.bw meth/CHG.bw meth/CHH.bw
exit 0
