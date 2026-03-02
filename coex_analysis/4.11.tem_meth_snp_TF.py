#!/usr/bin/env python3

"""
分析intersection.csv文件中的位点与TF motif区域的关系
"""

import pandas as pd
import os
from pathlib import Path

# 设置基础路径
BASE_DIR = Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana")
ANALYSIS_DIR = BASE_DIR / "list/xgboot/TPM_4.5/feature_analysis"
FEATURE_DIR = BASE_DIR / "list/xgboot/TPM_4.5/feature_data_extraction"

# TF文件列表（正样本和负样本）
TF_FILES = [
    BASE_DIR / "list/coex2/total/gene.TF_ids.txt",
    BASE_DIR / "list/coex2/10C/gene.TF_ids.txt",
    BASE_DIR / "list/coex2/16C/gene.TF_ids.txt",
    BASE_DIR / "list/coex2/22C/gene.TF_ids.txt",
    BASE_DIR / "list/coex2/total/gene.TF_negative_ids.txt",
    BASE_DIR / "list/coex2/10C/gene.TF_negative_ids.txt",
    BASE_DIR / "list/coex2/16C/gene.TF_negative_ids.txt",
    BASE_DIR / "list/coex2/22C/gene.TF_negative_ids.txt",
]

# 染色体名称映射 (NC_xxx -> Chr格式)
CHR_MAPPING = {
    'NC_003070.9': 'Chr1',
    'NC_003071.7': 'Chr2',
    'NC_003074.8': 'Chr3',
    'NC_003075.7': 'Chr4',
    'NC_003076.8': 'Chr5',
}

# 所有可能的feature_data_extraction子文件夹
FEATURE_FOLDERS = ["A_per_TPM", "A_TPM", "B_per_TPM", "B_TPM", "total_TPM"]

def read_all_tf_data(tf_files):
    """
    读取多个TF motif数据文件并合并
    返回包含motif_id, chr, start, stop的DataFrame
    """
    all_dfs = []
    
    for tf_file in tf_files:
        if not tf_file.exists():
            print(f"警告: TF文件不存在 {tf_file}")
            continue
        
        print(f"读取TF文件: {tf_file}")
        try:
            df = pd.read_csv(tf_file, sep='\t')
            # 选择需要的列并重命名
            df = df[['motif_id', 'sequence_name', 'start', 'stop']].copy()
            df.rename(columns={'sequence_name': 'chr'}, inplace=True)
            
            # 转换染色体名称
            df['chr'] = df['chr'].map(CHR_MAPPING)
            # 删除未映射的染色体
            df = df.dropna(subset=['chr'])
            
            print(f"  读取了 {len(df)} 个TF motif记录")
            all_dfs.append(df)
        except Exception as e:
            print(f"  错误: 读取文件失败 {tf_file.name}: {e}")
    
    if not all_dfs:
        print("错误: 没有成功读取任何TF文件！")
        return pd.DataFrame(columns=['motif_id', 'chr', 'start', 'stop'])
    
    # 合并所有数据并去重
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df = combined_df.drop_duplicates()
    
    print(f"\n总计读取了 {len(combined_df)} 个TF motif记录（去重后）")
    print(f"染色体分布: {combined_df['chr'].value_counts().to_dict()}")
    return combined_df

def load_all_feature_positions():
    """
    加载所有feature_data_extraction文件夹中的位置信息
    返回dict: {Site_ID: {'chr': chr, 'start': start, 'end': end}}
    """
    print("\n加载所有特征位置数据...")
    all_positions = {}
    
    for folder in FEATURE_FOLDERS:
        folder_path = FEATURE_DIR / folder
        if not folder_path.exists():
            print(f"  警告: 文件夹不存在 {folder_path}")
            continue
        
        print(f"  读取文件夹: {folder}")
        
        # 读取每种类型的数据
        for feature_type in ['CG', 'CHG', 'CHH', 'SNP']:
            if feature_type == 'SNP':
                file_path = folder_path / f"{folder}_snp_data.csv"
            else:
                file_path = folder_path / f"{folder}_{feature_type}_data.csv"
            
            if not file_path.exists():
                continue
            
            try:
                df = pd.read_csv(file_path)
                
                if feature_type == 'SNP':
                    # SNP文件格式: Feature,#Chromosome,Position,...
                    if 'Feature' in df.columns and '#Chromosome' in df.columns and 'Position' in df.columns:
                        for _, row in df.iterrows():
                            site_id = row['Feature']
                            if site_id not in all_positions:
                                all_positions[site_id] = {
                                    'chr': row['#Chromosome'],
                                    'start': row['Position'],
                                    'end': row['Position']
                                }
                else:
                    # 甲基化文件格式: Feature,chr,start,end,...
                    if 'Feature' in df.columns and 'chr' in df.columns and 'start' in df.columns:
                        for _, row in df.iterrows():
                            site_id = row['Feature']
                            if site_id not in all_positions:
                                # 甲基化使用start+1
                                all_positions[site_id] = {
                                    'chr': row['chr'],
                                    'start': row['start'] + 1,
                                    'end': row['start'] + 1
                                }
            except Exception as e:
                print(f"    警告: 读取文件失败 {file_path.name}: {e}")
    
    print(f"  总计加载了 {len(all_positions)} 个位置记录")
    return all_positions

def find_overlapping_motifs(site_chr, site_pos, tf_df):
    """
    找到与位点重叠的所有TF motif
    
    Parameters:
    -----------
    site_chr: str, 染色体名称
    site_pos: int, 位点位置
    tf_df: DataFrame, TF motif数据
    
    Returns:
    --------
    list of dict, 每个dict包含motif_id, motif_start, motif_stop
    """
    # 筛选同一染色体的motif
    chr_motifs = tf_df[tf_df['chr'] == site_chr]
    
    # 找到重叠的motif (site_pos在motif的start往前10bp和stop往后10bp之间)
    overlapping = chr_motifs[
        (chr_motifs['start'] - 10 <= site_pos) & 
        (chr_motifs['stop'] + 10 >= site_pos)
    ]
    
    if len(overlapping) == 0:
        return []
    
    # 返回所有重叠的motif信息
    results = []
    for _, row in overlapping.iterrows():
        results.append({
            'motif_id': row['motif_id'],
            'motif_start': row['start'],
            'motif_stop': row['stop']
        })
    
    return results

def process_intersection_file(intersection_file, all_positions, tf_df):
    """
    处理单个intersection.csv文件，添加位置和TF信息
    
    Parameters:
    -----------
    intersection_file: Path, intersection.csv文件路径
    all_positions: dict, 所有Site_ID的位置信息
    tf_df: DataFrame, TF motif数据
    
    Returns:
    --------
    DataFrame, 添加了位置和TF信息的结果
    """
    # 读取intersection.csv文件（只有Site_ID列表）
    df = pd.read_csv(intersection_file, header=None, names=['Site_ID'])
    # 删除第一行如果它是标题（不包含特征ID格式：CG_、CHG_、CHH_、SNP_）
    first_id = str(df.iloc[0]['Site_ID'])
    if not (first_id.startswith('CG_') or first_id.startswith('CHG_') or 
            first_id.startswith('CHH_') or first_id.startswith('SNP_')):
        df = df.iloc[1:].reset_index(drop=True)
    
    print(f"  原始记录数: {len(df)}")
    
    # 添加位置信息
    chr_list = []
    start_list = []
    end_list = []
    
    for site_id in df['Site_ID']:
        if site_id in all_positions:
            pos = all_positions[site_id]
            chr_list.append(pos['chr'])
            # 转换为整数
            start_list.append(int(pos['start']))
            end_list.append(int(pos['end']))
        else:
            chr_list.append(None)
            start_list.append(None)
            end_list.append(None)
    
    df['chr'] = chr_list
    # 使用 Int64 类型以支持缺失值的整数列
    df['position'] = pd.array(start_list, dtype='Int64')
    
    # 统计找到位置信息的数量
    n_with_pos = df['chr'].notna().sum()
    print(f"  找到位置信息: {n_with_pos}/{len(df)}")
    
    # 为每个位点找到重叠的TF motif
    all_overlapping_motifs = []
    
    for _, row in df.iterrows():
        if pd.notna(row['chr']) and pd.notna(row['position']):
            overlapping_motifs = find_overlapping_motifs(
                row['chr'], row['position'], tf_df
            )
            all_overlapping_motifs.append(overlapping_motifs)
        else:
            all_overlapping_motifs.append([])
    
    # 找出最大重叠motif数量
    max_motifs = max([len(motifs) for motifs in all_overlapping_motifs]) if all_overlapping_motifs else 0
    
    # 创建motif列的列名列表
    motif_cols = []
    
    if max_motifs > 0:
        print(f"  最大重叠motif数量: {max_motifs}")
        
        # 为每个motif创建单独的列
        for i in range(1, max_motifs + 1):
            motif_id_col = []
            motif_start_col = []
            motif_stop_col = []
            
            for motifs in all_overlapping_motifs:
                if i <= len(motifs):
                    # 有第i个motif
                    motif_id_col.append(motifs[i-1]['motif_id'])
                    motif_start_col.append(motifs[i-1]['motif_start'])
                    motif_stop_col.append(motifs[i-1]['motif_stop'])
                else:
                    # 没有第i个motif，填充NA
                    motif_id_col.append('NA')
                    motif_start_col.append('NA')
                    motif_stop_col.append('NA')
            
            df[f'motif_id_{i}'] = motif_id_col
            df[f'motif_start_{i}'] = motif_start_col
            df[f'motif_stop_{i}'] = motif_stop_col
            motif_cols.extend([f'motif_id_{i}', f'motif_start_{i}', f'motif_stop_{i}'])
    else:
        # 如果没有找到任何motif，添加一列NA
        df['motif_id_1'] = 'NA'
        df['motif_start_1'] = 'NA'
        df['motif_stop_1'] = 'NA'
        motif_cols = ['motif_id_1', 'motif_start_1', 'motif_stop_1']
    
    # 重新排列列顺序: Site_ID后面跟chr, position, 然后是其他列，最后是motif信息
    other_cols = [col for col in df.columns if col not in 
                  ['Site_ID', 'chr', 'position'] + motif_cols]
    new_order = ['Site_ID', 'chr', 'position'] + other_cols + motif_cols
    df = df[new_order]
    
    # 统计有TF motif的位点数量
    n_with_motif = sum(df['motif_id_1'] != 'NA')
    print(f"  位点在TF motif区域内: {n_with_motif}/{len(df)}")
    
    return df

def main():
    """
    主函数：处理feature_analysis文件夹下的intersection.csv文件
    """
    print("=" * 80)
    print("开始处理intersection.csv文件，添加位置和TF motif信息")
    print("=" * 80)
    
    # 读取所有TF motif数据
    tf_df = read_all_tf_data(TF_FILES)
    
    # 加载所有特征位置数据
    all_positions = load_all_feature_positions()
    
    # 定义要查找的文件后缀
    suffixes = ['_SA_intersection.csv', '_SB_intersection.csv', '_SPα_intersection.csv', 
                '_SPβ_intersection.csv', '_ST_intersection.csv']
    
    # 查找所有符合条件的文件
    all_files = []
    for suffix in suffixes:
        files = list(ANALYSIS_DIR.glob(f"*{suffix}"))
        all_files.extend(files)
    
    if not all_files:
        print(f"警告: 在 {ANALYSIS_DIR} 中未找到intersection.csv文件")
        return
    
    print(f"\n找到 {len(all_files)} 个intersection.csv文件")
    print("-" * 80)
    
    # 处理每个文件
    for intersection_file in sorted(all_files):
        filename = intersection_file.name
        print(f"\n处理文件: {filename}")
        
        # 处理文件
        result_df = process_intersection_file(intersection_file, all_positions, tf_df)
        
        # 保存结果 (覆盖原文件)
        output_file = intersection_file
        result_df.to_csv(output_file, index=False)
        print(f"  结果已保存到: {output_file.name}")
    
    print("\n" + "=" * 80)
    print("所有文件处理完成！")
    print("=" * 80)

if __name__ == "__main__":
    main()
