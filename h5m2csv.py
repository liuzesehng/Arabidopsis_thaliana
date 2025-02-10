#!/usr/bin/env python3

#
# h5m2csv.py: Convert HDF5 SNP matrix to CSV
#
# (c) 2021 by Joffrey Fitz (joffrey.fitz@tuebingen.mpg.de),
# Max Planck Institute for Developmental Biology,
# Tuebingen, Germany
#

import h5py
import numpy
import csv

f = h5py.File('/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/snp/1001_SNP_MATRIX/imputed_snps_binary.hdf5','r')

# Get all SNP positions for all chromosomes (len=10709949)
positions = f['positions'][:]

# Array of tupels with start/stop indices for each chromosome
chr_regions = f['positions'].attrs['chr_regions']

# Array of SNP positions for all chromosomes, each chromosome is a hash
# with "Chr<N>" as key, and a numpy.array of positions as value. 
snp_pos_on_chrs = [
	{ "label": "Chr1", "chr_idx": 0, "positions": positions[chr_regions[0][0]:chr_regions[0][1]] },
	{ "label": "Chr2", "chr_idx": 1, "positions": positions[chr_regions[1][0]:chr_regions[1][1]] },
	{ "label": "Chr3", "chr_idx": 2, "positions": positions[chr_regions[2][0]:chr_regions[2][1]] },
	{ "label": "Chr4", "chr_idx": 3, "positions": positions[chr_regions[3][0]:chr_regions[3][1]] },
	{ "label": "Chr5", "chr_idx": 4, "positions": positions[chr_regions[4][0]:chr_regions[4][1]] }
]

# 打开CSV文件进行写入
with open('snp_matrix.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)

    # 写入表头
    header = ["#Chromosome", "Position", "Count_zeros", "Count_ones"] + list(f['accessions'][:].astype(str))
    csvwriter.writerow(header)

    # 遍历所有染色体
    for chr in snp_pos_on_chrs:
        # 遍历所有位置信息
        for pos in numpy.nditer(chr["positions"]):
            # 找到特定位置的索引
            ix = numpy.where(chr["positions"] == pos)[0][0]

            # 将染色体起始位置加到SNP位置索引中
            ix = ix + chr_regions[chr["chr_idx"]][0]

            # 获取该位置对应的SNP数据
            snps = f['snps'][ix]

            # 统计0和1的数量
            cnt_zeros = numpy.count_nonzero(snps == 0)
            cnt_ones = numpy.count_nonzero(snps == 1)

            # 写入CSV文件
            row = [chr["label"], pos, cnt_zeros, cnt_ones] + list(snps.astype(str))
            csvwriter.writerow(row)

print("SNP矩阵已成功保存为'snp_matrix.csv'文件。")