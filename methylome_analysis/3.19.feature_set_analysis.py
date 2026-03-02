#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
特征集合分析脚本
从XGBoost优化结果中提取特征，进行集合运算并绘制韦恩图
同时分析RCA结果中的significant_sites.tsv文件，根据trend列进行分类分析
"""

import re
import pandas as pd
import os
import matplotlib.pyplot as plt
from upsetplot import plot as upset_plot
import glob

def extract_features_from_file(filepath):
    """从文件中提取特征列表"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 查找 "特征列表:" 后面的内容
        match = re.search(r'特征列表:\s*([^\n]+)', content)
        if match:
            feature_str = match.group(1)
            # 分割特征，去除空格
            features = [f.strip() for f in feature_str.split(',')]
            return set(features)
        else:
            print(f"警告: 在 {filepath} 中未找到特征列表")
            return set()
    except FileNotFoundError:
        print(f"错误: 文件不存在 - {filepath}")
        return set()
    except Exception as e:
        print(f"错误: 读取文件 {filepath} 时出错 - {str(e)}")
        return set()

def extract_sites_by_trend(rca_base_path, trend_category):
    """从RCA文件夹中提取特定trend类别的所有site"""
    sites = set()
    folders = ['A_TPM', 'A_per_TPM', 'B_TPM', 'B_per_TPM', 'total_TPM']
    
    for folder in folders:
        folder_path = os.path.join(rca_base_path, folder)
        if not os.path.exists(folder_path):
            print(f"  警告: 文件夹不存在 - {folder_path}")
            continue
        
        # 查找所有以significant_sites.tsv结尾的文件
        pattern = os.path.join(folder_path, '*significant_sites.tsv')
        files = glob.glob(pattern)
        
        for file in files:
            try:
                df = pd.read_csv(file, sep='\t')
                if 'trend' in df.columns and 'site' in df.columns:
                    # 提取特定trend的site
                    trend_sites = df[df['trend'] == trend_category]['site'].tolist()
                    sites.update(trend_sites)
            except Exception as e:
                print(f"  警告: 读取文件 {file} 时出错 - {str(e)}")
    
    return sites

def save_features_to_csv(feature_set, output_file, column_name='Feature'):
    """将特征集合保存为CSV文件"""
    df = pd.DataFrame(sorted(list(feature_set)), columns=[column_name])
    df.to_csv(output_file, index=False)
    print(f"  保存CSV: {output_file} ({len(feature_set)} 个特征)")

def main():
    # 基础路径
    base_path = "../../list/xgboot/TPM_4.5"
    output_dir = "../../list/xgboot/TPM_4.5/feature_analysis"
    rca_base_path = "../../list/RCA"
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print("特征集合分析 (XGBoost + RCA)")
    print("=" * 80)
    
    # ========================================
    # 1. 提取特征
    # ========================================
    print("\n[步骤 1] 提取特征列表...")
    
    # 表达量相关特征
    SA = extract_features_from_file(f"{base_path}/A_TPM_optimization_results.txt")
    SB = extract_features_from_file(f"{base_path}/B_TPM_optimization_results.txt")
    ST = extract_features_from_file(f"{base_path}/total_TPM_optimization_results.txt")
    
    print(f"  SA (A_TPM): {len(SA)} 个特征")
    print(f"  SB (B_TPM): {len(SB)} 个特征")
    print(f"  ST (total_TPM): {len(ST)} 个特征")
    
    # 剪接相关特征
    SPa = extract_features_from_file(f"{base_path}/A_per_TPM_optimization_results.txt")
    SPb = extract_features_from_file(f"{base_path}/B_per_TPM_optimization_results.txt")
    
    print(f"  SPα (A_per_TPM): {len(SPa)} 个特征")
    print(f"  SPβ (B_per_TPM): {len(SPb)} 个特征")
    
    # ========================================
    # 2. 计算核心集合
    # ========================================
    print("\n[步骤 2] 计算核心集合...")
    
    # Expression = SA ∩ SB ∩ ST (表达量核心)
    Expression = SA & SB & ST
    print(f"  Expression (SA ∩ SB ∩ ST): {len(Expression)} 个特征")
    
    # Expression_Ratio = SPα ∩ SPβ (表达比例核心)
    Expression_Ratio = SPa & SPb
    print(f"  Expression_Ratio (SPα ∩ SPβ): {len(Expression_Ratio)} 个特征")
    
    # ========================================
    # 3. 保存CSV文件
    # ========================================
    print("\n[步骤 3] 保存CSV文件...")
    
    # 保存核心集合
    save_features_to_csv(Expression, 
                        f"{output_dir}/Expression.csv",
                        "Expression_Feature")
    
    save_features_to_csv(Expression_Ratio,
                        f"{output_dir}/Expression_Ratio.csv",
                        "Expression_Ratio_Feature")
    
    # ========================================
    # 4. 生成Upset图
    # ========================================
    print("\n[步骤 4] 生成Upset图...")
    
    # 准备Upset图数据 - 使用5个原始集合
    all_features = SA | SB | ST | SPa | SPb
    
    # 创建布尔型DataFrame
    upset_data = pd.DataFrame({
        'ST': [f in ST for f in all_features],
        'SA': [f in SA for f in all_features],
        'SB': [f in SB for f in all_features],
        'SPα': [f in SPa for f in all_features],
        'SPβ': [f in SPb for f in all_features]
    }, index=list(all_features))
    
    # 转换为MultiIndex格式
    upset_data = upset_data.set_index(['ST', 'SA', 'SB', 'SPα', 'SPβ'])
    upset_counts = upset_data.groupby(level=[0, 1, 2, 3, 4]).size()
    
    # 绘制Upset图
    fig = plt.figure(figsize=(12, 7))
    upset_plot(upset_counts, 
               fig=fig,
               show_counts=True,
               element_size=40)
    
    plt.suptitle('Feature Set Intersection Analysis (ST, SA, SB, SPα, SPβ)', 
                 fontsize=14, 
                 fontweight='bold',
                 y=0.98)
    
    # 保存图片
    upset_file = f"{output_dir}/feature_upset_plot.png"
    plt.tight_layout()
    plt.savefig(upset_file, dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/feature_upset_plot.pdf", bbox_inches='tight')
    print(f"  保存Upset图: {upset_file}")
    print(f"  保存Upset图: {output_dir}/feature_upset_plot.pdf")
    plt.close()
    
    # ========================================
    # 5. 分析RCA结果 - 根据trend分类
    # ========================================
    print("\n[步骤 5] 分析RCA结果 - 提取trend分类的sites...")
    
    # 定义6个trend类别
    trend_categories = ['Up', 'Down', 'Up-like_HighTemp', 'Up-like_LowMid', 
                       'Down-like_HighTemp', 'Down-like_LowTemp']
    
    # 提取每个trend类别的sites
    trend_sets = {}
    for trend in trend_categories:
        sites = extract_sites_by_trend(rca_base_path, trend)
        trend_sets[trend] = sites
        print(f"  {trend}: {len(sites)} 个sites")
    
    # ========================================
    # 6. 生成Trend数据集与XGBoost特征集的Upset图
    # ========================================
    print("\n[步骤 6] 生成Trend与XGBoost特征集的Upset图...")
    
    # 合并所有集合（6个trend + 5个XGBoost特征集）
    all_items = (trend_sets['Up'] | trend_sets['Down'] | 
                 trend_sets['Up-like_HighTemp'] | trend_sets['Up-like_LowMid'] |
                 trend_sets['Down-like_HighTemp'] | trend_sets['Down-like_LowTemp'] |
                 ST | SA | SB | SPa | SPb)
    
    # 创建布尔型DataFrame
    upset_data_combined = pd.DataFrame({
        'Up': [item in trend_sets['Up'] for item in all_items],
        'Down': [item in trend_sets['Down'] for item in all_items],
        'Up-like_HighTemp': [item in trend_sets['Up-like_HighTemp'] for item in all_items],
        'Up-like_LowMid': [item in trend_sets['Up-like_LowMid'] for item in all_items],
        'Down-like_HighTemp': [item in trend_sets['Down-like_HighTemp'] for item in all_items],
        'Down-like_LowTemp': [item in trend_sets['Down-like_LowTemp'] for item in all_items],
        'ST': [item in ST for item in all_items],
        'SA': [item in SA for item in all_items],
        'SB': [item in SB for item in all_items],
        'SPα': [item in SPa for item in all_items],
        'SPβ': [item in SPb for item in all_items]
    }, index=list(all_items))
    
    # 转换为MultiIndex格式
    upset_data_combined = upset_data_combined.set_index([
        'Up', 'Down', 'Up-like_HighTemp', 'Up-like_LowMid', 
        'Down-like_HighTemp', 'Down-like_LowTemp',
        'ST', 'SA', 'SB', 'SPα', 'SPβ'
    ])
    upset_counts_combined = upset_data_combined.groupby(level=list(range(11))).size()
    
    # 绘制Upset图
    fig = plt.figure(figsize=(16, 9))
    upset_plot(upset_counts_combined, 
               fig=fig,
               show_counts=True,
               element_size=35)
    
    plt.suptitle('Trend Categories vs XGBoost Features Intersection Analysis', 
                 fontsize=14, 
                 fontweight='bold',
                 y=0.98)
    
    # 保存图片
    upset_combined_file = f"{output_dir}/trend_xgboost_upset_plot.png"
    plt.tight_layout()
    plt.savefig(upset_combined_file, dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/trend_xgboost_upset_plot.pdf", bbox_inches='tight')
    print(f"  保存Upset图: {upset_combined_file}")
    print(f"  保存Upset图: {output_dir}/trend_xgboost_upset_plot.pdf")
    plt.close()
    
    # ========================================
    # 7. 计算Trend类别与XGBoost特征集的交集
    # ========================================
    print("\n[步骤 7] 计算Trend类别与ST, SA, SB, SPα, SPβ的交集...")
    
    # 定义XGBoost特征集
    xgboost_sets = {
        'ST': ST,
        'SA': SA,
        'SB': SB,
        'SPα': SPa,
        'SPβ': SPb
    }
    
    # 存储所有交集结果
    trend_xgboost_intersections = {}
    
    for trend in trend_categories:
        print(f"  {trend}:")
        trend_xgboost_intersections[trend] = {}
        
        for xgb_name, xgb_set in xgboost_sets.items():
            # 计算交集
            intersection = trend_sets[trend] & xgb_set
            trend_xgboost_intersections[trend][xgb_name] = intersection
            
            print(f"    与{xgb_name}交集: {len(intersection)} 个特征")
            
            # 保存交集CSV文件（只保存非空交集）
            if len(intersection) > 0:
                save_features_to_csv(
                    intersection,
                    f"{output_dir}/{trend}_{xgb_name}_intersection.csv",
                    f"{trend}_{xgb_name}"
                )
    
    # ========================================
    # 8. 生成完整统计报告
    # ========================================
    print("\n[步骤 8] 生成完整统计报告...")
    
    # ========================================
    # 8. 生成完整统计报告
    # ========================================
    print("\n[步骤 8] 生成完整统计报告...")
    
    # 构建trend交集报告部分
    trend_report = "\n5. Trend类别统计:\n"
    for trend in trend_categories:
        trend_report += f"   - {trend}: {len(trend_sets[trend])} 个sites\n"
    
    trend_report += "\n6. Trend类别与XGBoost特征集的交集:\n"
    for trend in trend_categories:
        trend_report += f"   {trend}:\n"
        for xgb_name in ['ST', 'SA', 'SB', 'SPα', 'SPβ']:
            intersection_count = len(trend_xgboost_intersections[trend][xgb_name])
            trend_report += f"     - 与{xgb_name}交集: {intersection_count} 个特征\n"
    
    report = f"""
特征集合分析报告 (XGBoost + RCA)
{'=' * 80}

1. XGBoost原始特征集合:
   - SA (A_TPM): {len(SA)} 个特征
   - SB (B_TPM): {len(SB)} 个特征
   - ST (total_TPM): {len(ST)} 个特征
   - SPα (A_per_TPM): {len(SPa)} 个特征
   - SPβ (B_per_TPM): {len(SPb)} 个特征

2. XGBoost核心集合:
   - Expression (SA ∩ SB ∩ ST): {len(Expression)} 个特征
   - Expression_Ratio (SPα ∩ SPβ): {len(Expression_Ratio)} 个特征

3. XGBoost交集分析:
   - 总特征数: {len(all_features)} 个
   - ST独有: {len(ST - SA - SB - SPa - SPb)} 个特征
   - SA独有: {len(SA - ST - SB - SPa - SPb)} 个特征
   - SB独有: {len(SB - ST - SA - SPa - SPb)} 个特征
   - SPα独有: {len(SPa - ST - SA - SB - SPb)} 个特征
   - SPβ独有: {len(SPb - ST - SA - SB - SPa)} 个特征
   - Expression (SA ∩ SB ∩ ST): {len(Expression)} 个特征
   - Expression_Ratio (SPα ∩ SPβ): {len(Expression_Ratio)} 个特征
   - Expression与Expression_Ratio共同: {len(Expression & Expression_Ratio)} 个特征

4. RCA数据统计:
   - 总item数（Trend + XGBoost）: {len(all_items)} 个
{trend_report}

7. 输出文件:
   XGBoost核心集合CSV文件:
   - {output_dir}/Expression.csv
   - {output_dir}/Expression_Ratio.csv
   
   Trend与XGBoost特征集交集CSV文件:"""
    
    for trend in trend_categories:
        for xgb_name in ['ST', 'SA', 'SB', 'SPα', 'SPβ']:
            if len(trend_xgboost_intersections[trend][xgb_name]) > 0:
                report += f"\n   - {output_dir}/{trend}_{xgb_name}_intersection.csv"
    
    report += f"""
   
   可视化文件:
   - {output_dir}/feature_upset_plot.png (XGBoost特征集)
   - {output_dir}/feature_upset_plot.pdf
   - {output_dir}/trend_xgboost_upset_plot.png (Trend + XGBoost)
   - {output_dir}/trend_xgboost_upset_plot.pdf

{'=' * 80}
"""
    
    # 保存报告
    report_file = f"{output_dir}/analysis_report.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(report)
    print(f"\n报告已保存到: {report_file}")
    print("\n分析完成！")

if __name__ == "__main__":
    main()
