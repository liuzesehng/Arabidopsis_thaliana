#!/usr/bin/env python3

"""
分析summary_table文件中的位点与TF motif区域的关系
"""

import pandas as pd
import os
from pathlib import Path

# 设置基础路径
BASE_DIR = Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana")
RCA_DIR = BASE_DIR / "list/RCA"
FEATURE_DIR = BASE_DIR / "list/xgboot/TPM_4.5/feature_data_extraction"
TF_FILE = BASE_DIR / "list/coex2/total/gene.TF_ids.txt"

# 定义要处理的文件夹
FOLDERS = ["Scoupled_specific", "Sexpr_only", "Ssplice_only"]

# 染色体名称映射 (NC_xxx -> Chr格式)
CHR_MAPPING = {
    'NC_003071.7': 'Chr2',
}

def read_tf_data(tf_file):
    """
    读取TF motif数据文件
    返回包含motif_id, chr, start, stop的DataFrame
    """
    print(f"读取TF文件: {tf_file}")
    df = pd.read_csv(tf_file, sep='\t')
    # 选择需要的列并重命名
    df = df[['motif_id', 'sequence_name', 'start', 'stop']].copy()
    df.rename(columns={'sequence_name': 'chr'}, inplace=True)
    
    # 转换染色体名称
    df['chr'] = df['chr'].map(CHR_MAPPING)
    # 删除未映射的染色体
    df = df.dropna(subset=['chr'])
    
    print(f"读取了 {len(df)} 个TF motif记录")
    print(f"染色体分布: {df['chr'].value_counts().to_dict()}")
    return df

def read_feature_positions(folder_name, feature_type):
    """
    从feature_data_extraction文件夹读取位置信息
    
    Parameters:
    -----------
    folder_name: str, 文件夹名称
    feature_type: str, 特征类型 (CG, CHG, CHH, SNP)
    
    Returns:
    --------
    DataFrame with columns: Site_ID, chr, start, end
    """
    if feature_type.upper() == 'SNP':
        file_path = FEATURE_DIR / folder_name / f"{folder_name}_snp_data.csv"
    else:
        file_path = FEATURE_DIR / folder_name / f"{folder_name}_{feature_type}_data.csv"
    
    if not file_path.exists():
        print(f"警告: 文件不存在 {file_path}")
        return None
    
    print(f"  读取特征位置文件: {file_path.name}")
    df = pd.read_csv(file_path)
    
    if feature_type.upper() == 'SNP':
        # SNP文件格式: Feature,#Chromosome,Position,...
        df = df[['Feature', '#Chromosome', 'Position']].copy()
        df.rename(columns={
            'Feature': 'Site_ID',
            '#Chromosome': 'chr',
            'Position': 'start'
        }, inplace=True)
        # SNP使用start作为位置
        df['end'] = df['start']
    else:
        # 甲基化文件格式: Feature,chr,start,end,...
        df = df[['Feature', 'chr', 'start', 'end']].copy()
        df.rename(columns={'Feature': 'Site_ID'}, inplace=True)
        # 甲基化使用start+1
        df['start'] = df['start'] + 1
        df['end'] = df['start']  # 单碱基位点
    
    return df

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
    
    # 找到重叠的motif (site_pos在motif的start和stop之间)
    overlapping = chr_motifs[
        (chr_motifs['start'] <= site_pos) & 
        (chr_motifs['stop'] >= site_pos)
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

def process_summary_table(summary_file, feature_positions, tf_df):
    """
    处理单个summary_table文件，添加位置和TF信息
    
    Parameters:
    -----------
    summary_file: Path, summary_table文件路径
    feature_positions: DataFrame, 特征位置信息
    tf_df: DataFrame, TF motif数据
    
    Returns:
    --------
    DataFrame, 添加了位置和TF信息的结果
    """
    print(f"  处理文件: {summary_file.name}")
    
    # 读取summary_table
    df = pd.read_csv(summary_file, sep='\t')
    print(f"    原始记录数: {len(df)}")
    
    # 删除已存在的位置和motif列（如果存在）
    cols_to_remove = ['chr', 'start', 'end']
    # 删除所有motif相关列（包括motif_id_1, motif_id_2等）
    motif_pattern_cols = [col for col in df.columns if col.startswith('motif_')]
    cols_to_remove.extend(motif_pattern_cols)
    df = df.drop(columns=[col for col in cols_to_remove if col in df.columns])
    
    # 合并位置信息
    df = df.merge(feature_positions, on='Site_ID', how='left')
    
    # 为每个位点找到重叠的TF motif
    all_overlapping_motifs = []
    
    for _, row in df.iterrows():
        if pd.notna(row['chr']) and pd.notna(row['start']):
            overlapping_motifs = find_overlapping_motifs(
                row['chr'], row['start'], tf_df
            )
            all_overlapping_motifs.append(overlapping_motifs)
        else:
            all_overlapping_motifs.append([])
    
    # 找出最大重叠motif数量
    max_motifs = max([len(motifs) for motifs in all_overlapping_motifs]) if all_overlapping_motifs else 0
    
    # 创建motif列的列名列表
    motif_cols = []
    
    if max_motifs > 0:
        print(f"    最大重叠motif数量: {max_motifs}")
        
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
    
    # 重新排列列顺序: Site_ID后面跟chr, start, end, 然后是其他列，最后是motif信息
    other_cols = [col for col in df.columns if col not in 
                  ['Site_ID', 'chr', 'start', 'end'] + motif_cols]
    new_order = ['Site_ID', 'chr', 'start', 'end'] + other_cols + motif_cols
    df = df[new_order]
    
    # 统计有TF motif的位点数量
    n_with_motif = sum(df['motif_id_1'] != 'NA')
    print(f"    位点在TF motif区域内: {n_with_motif}/{len(df)}")
    
    return df

def main():
    """
    主函数：处理所有文件夹和summary_table文件
    """
    print("=" * 80)
    print("开始处理summary_table文件，添加位置和TF motif信息")
    print("=" * 80)
    
    # 读取TF motif数据
    tf_df = read_tf_data(TF_FILE)
    
    # 处理每个文件夹
    for folder in FOLDERS:
        print(f"\n处理文件夹: {folder}")
        print("-" * 80)
        
        rca_folder = RCA_DIR / folder
        
        # 查找所有summary_table文件
        summary_files = list(rca_folder.glob("*_summary_table.tsv"))
        
        if not summary_files:
            print(f"  警告: 在 {rca_folder} 中未找到summary_table文件")
            continue
        
        print(f"  找到 {len(summary_files)} 个summary_table文件")
        
        # 处理每个summary_table文件
        for summary_file in summary_files:
            # 确定特征类型 (从文件名提取)
            filename = summary_file.name
            if filename.startswith('SNP_'):
                feature_type = 'SNP'
            elif filename.startswith('CG_'):
                feature_type = 'CG'
            elif filename.startswith('CHG_'):
                feature_type = 'CHG'
            elif filename.startswith('CHH_'):
                feature_type = 'CHH'
            else:
                print(f"  警告: 无法识别文件类型 {filename}")
                continue
            
            # 读取特征位置信息
            feature_positions = read_feature_positions(folder, feature_type)
            
            if feature_positions is None:
                print(f"  跳过文件 {filename} (无法读取位置信息)")
                continue
            
            # 处理summary_table文件
            result_df = process_summary_table(summary_file, feature_positions, tf_df)
            
            # 保存结果 (覆盖原文件)
            output_file = summary_file
            result_df.to_csv(output_file, sep='\t', index=False)
            print(f"    结果已保存到: {output_file}")
    
    print("\n" + "=" * 80)
    print("所有文件处理完成！")
    print("=" * 80)

if __name__ == "__main__":
    main()
