#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UpSet plot for 17 sets:
  - 6 expression categories (Up, Down, Up-like_HighTemp, Up-like_LowMid, Down-like_HighTemp, Down-like_LowTemp)
  - 5 dataset types (SA, SB, SPα, SPβ, ST) from lun/TPM_4.5 optimization results
  - 6 TF genes (target_genes)
17个集合的UpSet图：6个表达类别 + 5个数据集类型 + 6个TF基因
"""

import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
import os
from itertools import combinations

# 设置路径
base_dir = "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng"
feature_dir = os.path.join(base_dir, "Arabidopsis_thaliana/list/xgboot/TPM_4.5/feature_analysis")
output_dir = os.path.join(base_dir, "Arabidopsis_thaliana/list/xgboot/TPM_4.5/feature_data_TF_annotation")
# lun/TPM_4.5 目录（SA/SB/SPα/SPβ/ST 数据来源）
lun_base_path = os.path.join(base_dir, "Arabidopsis_thaliana/lun/TPM_4.5")

# 目标TF基因（作为6个独立的集合）
target_genes = ['AT2G22540', 'AT5G10140', 'AT1G01060', 'AT3G09600', 'AT5G02840', 'AT5G17300']

# 数据集类型（作为5个独立的集合）——数据来源改为 lun/TPM_4.5 优化结果文件
dataset_types = ['SA', 'SB', 'SPα', 'SPβ', 'ST']

# dataset_type → 对应的优化结果文件名
DATASET_FILE_MAP = {
    'SA':  'A_TPM_optimization_results.txt',
    'SB':  'B_TPM_optimization_results.txt',
    'ST':  'total_TPM_optimization_results.txt',
    'SPα': 'A_per_TPM_optimization_results.txt',
    'SPβ': 'B_per_TPM_optimization_results.txt',
}

# 表达类别（作为6个独立的集合）——根据 feature_analysis 目录中文件名前缀区分
expression_categories = ['Up', 'Down', 'Up-like_HighTemp', 'Up-like_LowMid',
                         'Down-like_HighTemp', 'Down-like_LowTemp']

def extract_features_from_file(file_path):
    """从 lun/TPM_4.5 优化结果 txt 文件中提取特征（Site_ID）列表"""
    site_ids = set()
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                if '特征列表:' in line:
                    parts = line.split('特征列表:', 1)
                    if len(parts) > 1:
                        features_str = parts[1].strip()
                        features = [feat.strip() for feat in features_str.split(',') if feat.strip()]
                        site_ids.update(features)
        print(f"  {os.path.basename(file_path)}: {len(site_ids)} 个位点")
    except Exception as e:
        print(f"  {file_path}: 读取出错 - {e}")
    return site_ids


def collect_dataset_type_sites(dataset_type):
    """从 lun/TPM_4.5 优化结果文件中收集某个数据集类型的 Site_ID"""
    file_name = DATASET_FILE_MAP.get(dataset_type)
    if not file_name:
        print(f"  {dataset_type}: 未找到对应文件映射")
        return set()
    file_path = os.path.join(lun_base_path, file_name)
    site_ids = extract_features_from_file(file_path)
    print(f"  {dataset_type} 总计: {len(site_ids)} 个位点")
    return site_ids


def collect_expression_category_sites(category):
    """收集某个表达类别（如 Up、Down-like_HighTemp 等）的所有 Site_ID
    
    匹配规则：文件名以 '{category}_' 开头，以 '_intersection.csv' 结尾，
    且不含 'target_genes' 或 'Expression' 关键字。
    """
    all_files = os.listdir(feature_dir)
    matching_files = [
        f for f in all_files
        if f.startswith(f'{category}_')
        and f.endswith('_intersection.csv')
        and 'target_genes' not in f
        and 'Expression' not in f
    ]
    
    all_site_ids = set()
    
    for file_name in matching_files:
        file_path = os.path.join(feature_dir, file_name)
        try:
            df = pd.read_csv(file_path)
            if not df.empty and 'Site_ID' in df.columns:
                site_ids = set(df['Site_ID'].tolist())
                all_site_ids.update(site_ids)
                print(f"  {file_name}: {len(site_ids)} 个位点")
        except Exception as e:
            print(f"  {file_name}: 读取出错 - {e}")
    
    print(f"  {category} 总计: {len(all_site_ids)} 个位点")
    return all_site_ids

def collect_gene_sites(gene_id):
    """收集包含某个基因的所有Site_ID"""
    all_files = os.listdir(feature_dir)
    # 只处理以这些后缀结尾的文件
    suffixes = ['_SA_intersection.csv', '_SB_intersection.csv', 
                '_SPα_intersection.csv', '_SPβ_intersection.csv', 
                '_ST_intersection.csv']
    
    all_site_ids = set()
    
    for file_name in all_files:
        # 检查是否是目标文件
        if not any(file_name.endswith(suffix) for suffix in suffixes):
            continue
            
        file_path = os.path.join(feature_dir, file_name)
        try:
            df = pd.read_csv(file_path)
            if df.empty:
                continue
            
            # 获取所有motif_id列
            motif_cols = [col for col in df.columns if 'motif_id' in col]
            
            # 检查每行是否包含目标基因
            for idx, row in df.iterrows():
                motif_values = [str(row[col]) for col in motif_cols]
                if gene_id in motif_values:
                    all_site_ids.add(row['Site_ID'])
            
        except Exception as e:
            print(f"  {file_name}: 处理出错 - {e}")
    
    print(f"  {gene_id}: {len(all_site_ids)} 个位点")
    return all_site_ids

def main():
    """主函数"""
    print("=" * 80)
    print("收集17个集合的数据（6个表达类别 + 5个数据集类型 + 6个TF基因）")
    print("=" * 80)
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 收集所有集合的数据
    all_sets = {}
    
    # 1. 收集6个表达类别的Site_ID（来自 feature_analysis 目录 CSV 文件前缀）
    print("\n1. 收集表达类别的Site_ID:")
    print("-" * 80)
    for category in expression_categories:
        print(f"\n处理 {category}:")
        site_ids = collect_expression_category_sites(category)
        all_sets[category] = site_ids
    
    # 2. 收集5个数据集类型的Site_ID（来自 lun/TPM_4.5 优化结果文件）
    print("\n2. 收集数据集类型的Site_ID:")
    print("-" * 80)
    for dataset_type in dataset_types:
        print(f"\n处理 {dataset_type}:")
        site_ids = collect_dataset_type_sites(dataset_type)
        all_sets[dataset_type] = site_ids
    
    # 3. 收集6个TF基因的Site_ID
    print("\n3. 收集TF基因的Site_ID:")
    print("-" * 80)
    for gene in target_genes:
        print(f"\n处理 {gene}:")
        site_ids = collect_gene_sites(gene)
        all_sets[gene] = site_ids
    
    # 4. 保存每个集合的Site_ID到文件
    print("\n4. 保存各集合的Site_ID:")
    print("-" * 80)
    for set_name, site_ids in all_sets.items():
        if site_ids:
            safe_name = set_name.replace('/', '_').replace(' ', '_')
            output_file = os.path.join(output_dir, f'{safe_name}_site_ids.txt')
            with open(output_file, 'w') as f:
                for site_id in sorted(site_ids):
                    f.write(f"{site_id}\n")
            print(f"  {set_name}: {len(site_ids)} 个位点 -> {safe_name}_site_ids.txt")
    
    # 5. 创建UpSet图
    print("\n5. 创建UpSet图:")
    print("-" * 80)
    create_upset_plot(all_sets)
    
    # 6. 保存交集统计
    print("\n6. 保存交集统计:")
    print("-" * 80)
    save_intersection_stats(all_sets)
    
    # 7. 保存交集的Site_ID
    print("\n7. 保存各种交集的Site_ID:")
    print("-" * 80)
    save_intersection_site_ids(all_sets)
    
    print("\n" + "=" * 80)
    print("分析完成！")
    print(f"所有文件已保存到: {output_dir}")
    print(f"集合数量: {len(all_sets)}")
    print("  表达类别:", expression_categories)
    print("  数据集类型:", dataset_types)
    print("  TF基因:", target_genes)
    print("=" * 80)

def create_upset_plot(data_dict):
    """创建9个集合的UpSet图"""
    # 准备数据
    memberships = []
    all_site_ids = set().union(*data_dict.values())
    
    print(f"  总共有 {len(all_site_ids)} 个不同的Site_ID")
    
    for site_id in all_site_ids:
        membership = tuple(set_name for set_name, site_ids in data_dict.items() 
                          if site_id in site_ids)
        memberships.append(membership)
    
    # 创建UpSet数据
    upset_data = from_memberships(memberships)
    
    # 创建图形
    fig = plt.figure(figsize=(16, 10))
    upset = UpSet(upset_data, 
                  subset_size='count',
                  show_counts=True,
                  element_size=50,
                  intersection_plot_elements=10)
    upset.plot(fig=fig)
    
    n_sets = len(data_dict)
    plt.suptitle(
        f'UpSet Plot: {n_sets} Sets\n'
        '6 Expression Categories (Up/Down/Up-like_HighTemp/Up-like_LowMid/Down-like_HighTemp/Down-like_LowTemp)\n'
        '+ 5 Dataset Types (SA/SB/SPα/SPβ/ST) + 6 TF Genes',
        fontsize=13, y=0.98
    )
    
    # 保存图形
    output_file = os.path.join(output_dir, 'upset_17_sets.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  PNG图已保存: upset_17_sets.png")
    
    output_file_pdf = os.path.join(output_dir, 'upset_17_sets.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight')
    print(f"  PDF图已保存: upset_17_sets.pdf")
    
    plt.close()

def save_intersection_stats(data_dict):
    """保存交集统计信息"""
    stats = []
    set_names = list(data_dict.keys())
    
    # 单独集合统计
    for name in set_names:
        stats.append({
            'Sets': name,
            'Count': len(data_dict[name]),
            'Type': 'Single',
            'Set_Size': 1
        })
    
    # 统计2到9个集合的所有交集
    print("\n  计算各级别交集...")
    intersection_summary = {}  # 用于统计每个级别的交集数量和元素数量
    
    n_total = len(set_names)
    for r in range(2, n_total + 1):  # 从2个集合到全部集合
        print(f"  计算 {r} 个集合的交集...")
        intersection_summary[r] = {'count': 0, 'total_elements': 0, 'details': []}  # type: ignore[index]
        
        for combo in combinations(set_names, r):
            # 计算交集
            intersection = set.intersection(*[data_dict[name] for name in combo])
            intersection_size = len(intersection)
            
            # 记录到stats
            sets_str = ' ∩ '.join(combo)
            stats.append({
                'Sets': sets_str,
                'Count': intersection_size,
                'Type': f'{r}-way',
                'Set_Size': r
            })
            
            # 统计汇总
            if intersection_size > 0:
                intersection_summary[r]['count'] += 1
                intersection_summary[r]['total_elements'] += intersection_size
                intersection_summary[r]['details'].append({
                    'sets': combo,
                    'size': intersection_size,
                    'site_ids': intersection
                })
    
    # 保存统计
    df_stats = pd.DataFrame(stats)
    output_file = os.path.join(output_dir, 'intersection_stats.csv')
    df_stats.to_csv(output_file, index=False)
    print(f"  交集统计已保存: intersection_stats.csv")
    
    # 保存汇总统计到单独文件
    summary_data = []
    for r in range(2, n_total + 1):
        summary_data.append({
            'Intersection_Level': f'{r} sets',
            'Number_of_Intersections': intersection_summary[r]['count'],
            'Total_Elements': intersection_summary[r]['total_elements'],
            'Non_Empty_Intersections': len([d for d in intersection_summary[r]['details'] if d['size'] > 0])
        })
    
    df_summary = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, 'intersection_level_summary.csv')
    df_summary.to_csv(summary_file, index=False)
    print(f"  交集级别汇总已保存: intersection_level_summary.csv")
    
    # 显示重要统计
    print(f"\n  关键统计:")
    print(f"    单集合: {len([s for s in stats if s['Type'] == 'Single'])} 个")
    for r in range(2, n_total + 1):
        non_empty = len([d for d in intersection_summary[r]['details'] if d['size'] > 0])
        total_elements = intersection_summary[r]['total_elements']
        print(f"    {r}集合交集: 共{intersection_summary[r]['count']}个组合, "
              f"非空{non_empty}个, 共{total_elements}个元素")
    
    # 保存每个级别的详细交集信息到单独文件
    for r in range(2, n_total + 1):
        if intersection_summary[r]['details']:
            detail_file = os.path.join(output_dir, f'intersection_{r}_way_details.txt')
            with open(detail_file, 'w') as f:
                f.write(f"# {r}-way Intersections\n")
                f.write(f"# Total combinations: {intersection_summary[r]['count']}\n")
                f.write(f"# Non-empty intersections: {len([d for d in intersection_summary[r]['details'] if d['size'] > 0])}\n")
                f.write(f"# Total elements: {intersection_summary[r]['total_elements']}\n\n")
                
                for detail in intersection_summary[r]['details']:
                    if detail['size'] > 0:  # 只保存非空交集
                        f.write(f"\n{'='*80}\n")
                        f.write(f"Sets: {' ∩ '.join(detail['sets'])}\n")
                        f.write(f"Count: {detail['size']} sites\n")
                        f.write(f"Site IDs:\n")
                        for site_id in sorted(detail['site_ids']):
                            f.write(f"  {site_id}\n")
            
            print(f"  {r}集合交集详情已保存: intersection_{r}_way_details.txt")
    
    return intersection_summary

def save_intersection_site_ids(data_dict):
    """保存重要交集的Site_ID到文件"""
    n_total = len(data_dict)
    
    # 保存全部集合交集
    all_intersection = set.intersection(*data_dict.values())
    if all_intersection:
        output_file = os.path.join(output_dir, f'intersection_all_{n_total}_sets.txt')
        with open(output_file, 'w') as f:
            f.write(f"# All {n_total} sets intersection: {len(all_intersection)} sites\n")
            for site_id in sorted(all_intersection):
                f.write(f"{site_id}\n")
        print(f"  全集合交集 ({len(all_intersection)} 个位点) -> intersection_all_{n_total}_sets.txt")
    
    # 保存6个表达类别的交集
    valid_cats = [c for c in expression_categories if c in data_dict and data_dict[c]]
    if len(valid_cats) >= 2:
        cat_intersection = set.intersection(*[data_dict[c] for c in valid_cats])
        if cat_intersection:
            output_file = os.path.join(output_dir, 'intersection_expression_categories.txt')
            with open(output_file, 'w') as f:
                f.write(f"# {len(valid_cats)} expression categories intersection: {len(cat_intersection)} sites\n")
                f.write(f"# Categories: {', '.join(valid_cats)}\n")
                for site_id in sorted(cat_intersection):
                    f.write(f"{site_id}\n")
            print(f"  {len(valid_cats)}个表达类别交集 ({len(cat_intersection)} 个位点) -> intersection_expression_categories.txt")
    
    # 保存5个数据集类型的交集
    valid_dts = [dt for dt in dataset_types if dt in data_dict and data_dict[dt]]
    if len(valid_dts) >= 2:
        dataset_intersection = set.intersection(*[data_dict[dt] for dt in valid_dts])
        if dataset_intersection:
            output_file = os.path.join(output_dir, 'intersection_5_datasets.txt')
            with open(output_file, 'w') as f:
                f.write(f"# {len(valid_dts)} dataset types intersection: {len(dataset_intersection)} sites\n")
                for site_id in sorted(dataset_intersection):
                    f.write(f"{site_id}\n")
            print(f"  {len(valid_dts)}个数据集类型交集 ({len(dataset_intersection)} 个位点) -> intersection_5_datasets.txt")
    
    # 保存6个TF基因的交集
    valid_genes = [g for g in target_genes if g in data_dict and data_dict[g]]
    if len(valid_genes) >= 2:
        gene_intersection = set.intersection(*[data_dict[g] for g in valid_genes])
        if gene_intersection:
            output_file = os.path.join(output_dir, 'intersection_6_genes.txt')
            with open(output_file, 'w') as f:
                f.write(f"# {len(valid_genes)} TF genes intersection: {len(gene_intersection)} sites\n")
                for site_id in sorted(gene_intersection):
                    f.write(f"{site_id}\n")
            print(f"  {len(valid_genes)}个TF基因交集 ({len(gene_intersection)} 个位点) -> intersection_6_genes.txt")

if __name__ == "__main__":
    main()
