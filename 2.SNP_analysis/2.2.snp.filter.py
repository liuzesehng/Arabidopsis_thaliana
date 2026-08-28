import pandas as pd

# 加载数据
data = pd.read_csv('snp_matrix.csv', sep=',')
# 根据 Position 列筛选数据（16568746 <= Position <= 16575692）
data = data[(data['Position'] >= 16568746) & (data['Position'] <= 16575692)]

# 计算 MAF
# Count_zeros 是第3列，Count_ones 是第4列
data['total_alleles'] = data['Count_zeros'] + data['Count_ones']

# 计算次等位基因频率（取较小值的频率）
data['ref_freq'] = data['Count_zeros'] / data['total_alleles']
data['alt_freq'] = data['Count_ones'] / data['total_alleles']
data['MAF'] = data[['ref_freq', 'alt_freq']].min(axis=1)

# 创建结果 DataFrame，包含染色体、位置和频率信息
maf_results = data[['#Chromosome', 'Position', 'Count_zeros', 'Count_ones', 
                     'total_alleles', 'ref_freq', 'alt_freq', 'MAF']]

# 过滤掉 MAF 小于 0.01 的 SNP
filtered_maf_results = maf_results[maf_results['MAF'] >= 0.01]

# 同时过滤原始数据矩阵，只保留 MAF >= 0.01 的行
filtered_data = data[data['MAF'] >= 0.01]
# 保存完整的 MAF 统计结果
maf_results.to_csv("snp_frequencies.csv", sep='\t', index=False)

# 保存过滤后的 MAF 统计结果
filtered_maf_results.to_csv("filtered_snp_frequencies.csv", sep='\t', index=False)

# 保存过滤后的完整 SNP 矩阵（包含所有样本列）
# 只保留原始列，去掉临时计算列
filtered_snp_matrix = filtered_data.iloc[:, :data.columns.get_loc('total_alleles')]
filtered_snp_matrix.to_csv("filtered_snp_matrix.csv", sep='\t', index=False)

print(f"Total SNPs: {len(data)}")
print(f"SNPs after MAF >= 0.01 filtering: {len(filtered_data)}")
print(f"Filtered out: {len(data) - len(filtered_data)} SNPs")
