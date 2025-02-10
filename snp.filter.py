import pandas as pd

# 加载数据
data = pd.read_csv('snp_results.txt', sep='\t')

# 统计每个position的SNP频率
maf_results = (
    data.groupby('position')['SNP']
    .value_counts(normalize=True)  # 计算频率
    .unstack(fill_value=0)  # 转为表格
    .rename(columns={0: 'ref_freq', 1: 'alt_freq'})  # 重命名列
)

# 计算MAF
maf_results['MAF'] = maf_results[['ref_freq', 'alt_freq']].min(axis=1)

# 过滤掉MAF小于0.05的SNP
filtered_maf_results = maf_results[maf_results['MAF'] >= 0.05]

# 保存到 CSV 文件
maf_results.to_csv("snp_frequencies.csv", sep='\t', index=True)

# 保存过滤后的结果到 CSV 文件
filtered_maf_results.to_csv("filtered_snp_frequencies.csv", sep='\t', index=True)
