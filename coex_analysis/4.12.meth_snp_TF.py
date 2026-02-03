#!/usr/bin/env python3

"""
为feature_data_extraction文件夹下的CSV文件添加TF motif位置信息
输出到新文件夹feature_data_TF_annotation
"""

import pandas as pd
import os
from pathlib import Path

# 设置基础路径
BASE_DIR = Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana")
FEATURE_DIR = BASE_DIR / "list/xgboot/TPM_4.5/feature_data_extraction"
OUTPUT_DIR = BASE_DIR / "list/xgboot/TPM_4.5/feature_data_TF_annotation"
TF_FILE = BASE_DIR / "list/coex2/total/gene.TF_ids.txt"

# 定义要处理的文件夹
FOLDERS = ["Scoupled_specific", "Sexpr_only", "Ssplice_only"]

# 染色体名称映射 (NC_xxx -> Chr格式)
# 拟南芥染色体映射
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
    # 删除未映射的染色体（如果有的话）
    df = df.dropna(subset=['chr'])
    
    print(f"读取了 {len(df)} 个TF motif记录")
    print(f"染色体分布: {df['chr'].value_counts().to_dict()}")
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

def process_csv_file(csv_file, file_type, tf_df, output_dir):
    """
    处理单个CSV文件，添加TF motif信息
    
    Parameters:
    -----------
    csv_file: Path, CSV文件路径
    file_type: str, 文件类型 ('snp', 'CG', 'CHG', 'CHH')
    tf_df: DataFrame, TF motif数据
    output_dir: Path, 输出目录
    
    Returns:
    --------
    DataFrame, 添加了TF信息的结果
    """
    print(f"  处理文件: {csv_file.name}")
    
    # 读取CSV文件
    df = pd.read_csv(csv_file)
    print(f"    原始记录数: {len(df)}")
    
    # 根据文件类型确定位点位置
    if file_type == 'snp':
        # SNP文件: 使用Position列
        df['site_chr'] = df['#Chromosome']
        df['site_pos'] = df['Position']
    else:
        # 甲基化文件: 使用start+1
        df['site_chr'] = df['chr']
        df['site_pos'] = df['start'] + 1
    
    # 为每个位点找到重叠的TF motif
    all_overlapping_motifs = []
    
    for _, row in df.iterrows():
        overlapping_motifs = find_overlapping_motifs(
            row['site_chr'], row['site_pos'], tf_df
        )
        all_overlapping_motifs.append(overlapping_motifs)
    
    # 找出最大重叠motif数量
    max_motifs = max([len(motifs) for motifs in all_overlapping_motifs]) if all_overlapping_motifs else 0
    
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
    else:
        # 如果没有找到任何motif，添加一列NA
        df['motif_id_1'] = 'NA'
        df['motif_start_1'] = 'NA'
        df['motif_stop_1'] = 'NA'
    
    # 删除临时列
    df = df.drop(columns=['site_chr', 'site_pos'])
    
    # 统计有TF motif的位点数量
    n_with_motif = sum(df['motif_id_1'] != 'NA')
    print(f"    位点在TF motif区域内: {n_with_motif}/{len(df)}")
    
    # 保存结果到新文件
    output_file = output_dir / csv_file.name
    df.to_csv(output_file, index=False)
    print(f"    结果已保存到: {output_file}")
    
    return df

def main():
    """
    主函数：处理所有CSV文件
    """
    print("=" * 80)
    print("开始处理feature CSV文件，添加TF motif信息")
    print("=" * 80)
    
    # 读取TF motif数据
    tf_df = read_tf_data(TF_FILE)
    
    # 创建输出根目录
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n输出目录: {OUTPUT_DIR}")
    
    # 处理每个文件夹
    for folder in FOLDERS:
        print(f"\n处理文件夹: {folder}")
        print("-" * 80)
        
        input_folder = FEATURE_DIR / folder
        output_folder = OUTPUT_DIR / folder
        
        # 创建输出文件夹
        output_folder.mkdir(parents=True, exist_ok=True)
        
        # 查找所有CSV文件
        csv_files = list(input_folder.glob("*.csv"))
        
        if not csv_files:
            print(f"  警告: 在 {input_folder} 中未找到CSV文件")
            continue
        
        print(f"  找到 {len(csv_files)} 个CSV文件")
        
        # 处理每个CSV文件
        for csv_file in csv_files:
            # 确定文件类型
            filename = csv_file.name
            if '_snp_data.csv' in filename:
                file_type = 'snp'
            elif '_CG_data.csv' in filename:
                file_type = 'CG'
            elif '_CHG_data.csv' in filename:
                file_type = 'CHG'
            elif '_CHH_data.csv' in filename:
                file_type = 'CHH'
            else:
                print(f"  警告: 无法识别文件类型 {filename}")
                continue
            
            # 处理CSV文件
            process_csv_file(csv_file, file_type, tf_df, output_folder)
    
    print("\n" + "=" * 80)
    print("所有文件处理完成！")
    print(f"结果保存在: {OUTPUT_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
