#!/bin/bash
source ~/.bashrc
conda activate bismark
for i in CG CHG CHH; do
    awk 'BEGIN{FS="\t"; OFS="\t"} {
        if ($4 == "") $4 = 0
        if ($5 == "" || NF < 5) $5 = 0
        print $1, $2, $3, $4, $5
    }' meth_experiment/${i}.CV.bedGraph > meth_experiment/${i}.CV.fill0.bedGraph
done
for i in CG CHG CHH
do
    cut -d $'\t' -f 1-3,5 meth_experiment/$i.CV.fill0.bedGraph > meth_experiment/$i.bedGraph
    bedGraphToBigWig meth_experiment/$i.bedGraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt meth_experiment/$i.bw
    rm -rf meth_experiment/$i.bedGraph
done
conda activate deeptools
computeMatrix scale-regions -S meth_experiment/CG.bw meth_experiment/CHG.bw meth_experiment/CHH.bw \
                            -R gene_RCA.bed \
                            --numberOfProcessors 10 \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --binSize 100 \
                            -o meth_experiment/meth.matrix.mat.gz

plotHeatmap -m meth_experiment/meth.matrix.mat.gz \
            --perGroup \
            --colorMap RdBu_r \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "RCA methylation"  \
            -out meth_experiment/meth.cv.pdf
rm -rf meth_experiment/CG.bw meth_experiment/CHG.bw meth_experiment/CHH.bw meth_experiment/meth.matrix.mat.gz

conda activate bismark
for i in CG CHG CHH
do
    cut -d $'\t' -f 1-3,4 meth_experiment/$i.CV.fill0.bedGraph > meth_experiment/$i.bedGraph
    bedGraphToBigWig meth_experiment/$i.bedGraph /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/faidx/Arabidopsis_thaliana.chrSize.txt meth_experiment/$i.bw
    rm -rf meth_experiment/$i.CV.fill0.bedGraph meth_experiment/$i.bedGraph 
done

conda activate deeptools

computeMatrix scale-regions -S meth_experiment/CG.bw meth_experiment/CHG.bw meth_experiment/CHH.bw \
                            -R gene_RCA.bed \
                            --numberOfProcessors 10 \
                            --beforeRegionStartLength 2000 \
                            --regionBodyLength 3000 \
                            --afterRegionStartLength 2000 \
                            --missingDataAsZero \
                            --binSize 100 \
                            -o meth_experiment/meth.matrix.mat.gz

plotHeatmap -m meth_experiment/meth.matrix.mat.gz \
            --perGroup \
            --colorMap RdBu_r \
            --heatmapHeight 7 \
            --heatmapWidth 15 \
            --regionsLabel "RCA methylation" \
            -out meth_experiment/meth.cv.pdf
rm -rf meth_experiment/CG.bw meth_experiment/CHG.bw meth_experiment/CHH.bw
