import pandas as pd
import os
import numpy as np

# 读取文件名文件
with open('filenames.txt', 'r') as f:
    filenames = [line.strip() for line in f]

# 循环处理每一个文件
for filename in filenames:
    # 读取文件
    df = pd.read_csv(filename, sep='\t')

    basename, ext = os.path.splitext(filename)
    new_filename = basename + "_count" + ext

    # 计算第四列的总和
    total_tpm = round(df['TPM'].sum(), 3)

    # 添加新的一列
    #df['SUM'] = 0
    df.loc[0, 'SUM'] = total_tpm

    # 计算第四列第二行与第三行的比值
    tpm_0 = df.loc[0, 'TPM']
    tpm_1 = df.loc[1, 'TPM']
    tpm_2 = df.loc[2, 'TPM']

    if pd.isna(tpm_0) or pd.isna(tpm_1) or tpm_1 == 0:
        ratio = np.nan
        print("Warning: invalid values encountered in TPM calculations.\n" + filename)
    else:
        ratio = round(tpm_0 / tpm_1, 3)
        ratio1 = round(tpm_0 / (tpm_1 + tpm_2), 3)
        ratio2 = round(tpm_1 / tpm_2, 3)
        level1 = round(tpm_0 / total_tpm, 3)
        level2 = round(tpm_1 / total_tpm, 3)
        level3 = round(tpm_2 / total_tpm, 3)


    # 添加新的一列
    #df['Ratio'] = 0
    df.loc[0, 'Ratio'] = ratio
    df.loc[1, 'Ratio'] = ratio1
    df.loc[2, 'Ratio'] = ratio2

    df.loc[0, 'Level'] = level1
    df.loc[1, 'Level'] = level2
    df.loc[2, 'Level'] = level3

    # 输出结果
    df.to_csv(f'{new_filename}', sep='\t', index=False)