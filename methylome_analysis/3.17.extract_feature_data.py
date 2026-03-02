#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
特征数据提取脚本
根据特征列表从SNP和甲基化数据文件中提取对应的行

特征格式说明:
- 特征名称格式: 类型_行号 (如 CHH_63, CG_43, snp_13)
- 甲基化文件 (CG/CHG/CHH): 无表头，行号从1开始
  例如: CHH_63 表示CHH文件的第63行
- SNP文件: 有表头，特征数字+1是实际行号
  例如: snp_13 表示SNP文件的第14行(包含表头)
"""

import pandas as pd
import os
import re

def read_feature_csv(filepath):
    """读取特征CSV文件，返回特征集合"""
    try:
        df = pd.read_csv(filepath)
        # 获取第一列的所有特征
        features = set(df.iloc[:, 0].values)
        print(f"  从 {os.path.basename(filepath)} 读取了 {len(features)} 个特征")
        return features
    except Exception as e:
        print(f"  错误: 读取 {filepath} 失败 - {str(e)}")
        return set()

def parse_feature_name(feature):
    """解析特征名称，提取类型和位置信息"""
    # 特征格式如: CHH_63, CG_43, snp_13, CHG_325
    match = re.match(r'([A-Za-z]+)_(\d+)', feature)
    if match:
        feature_type = match.group(1)
        position = int(match.group(2))
        return feature_type, position
    return None, None

def extract_snp_data(features, input_file, output_file):
    """从SNP文件中提取特征数据(基于行号,有表头)"""
    print(f"\n  处理SNP文件: {os.path.basename(input_file)}")
    
    # 筛选出snp相关的特征
    snp_features = {f for f in features if f.startswith('snp_')}
    if not snp_features:
        print(f"    没有找到snp特征")
        return
    
    print(f"    需要查找 {len(snp_features)} 个SNP特征")
    
    try:
        # 读取SNP文件（制表符分隔）
        df = pd.read_csv(input_file, sep='\t')
        
        # 提取行号（特征数字+1，因为有表头）
        rows_to_extract = []
        feature_names = []
        
        for feature in sorted(snp_features):
            _, pos = parse_feature_name(feature)
            if pos is not None:
                # SNP文件有表头，行号 = 特征数字 + 1，但在pandas中索引从0开始
                # 所以实际索引 = 特征数字
                if pos < len(df):
                    rows_to_extract.append(pos)
                    feature_names.append(feature)
        
        # 提取对应行
        if rows_to_extract:
            matched_df = df.iloc[rows_to_extract].copy()
            matched_df.insert(0, 'Feature', feature_names)
            
            # 按照start列（第3列）升序排序
            matched_df = matched_df.sort_values('Position')

            # 保存结果
            matched_df.to_csv(output_file, index=False)
            print(f"    -> 找到 {len(matched_df)} 行数据")
            print(f"    -> 保存到: {output_file}")
        else:
            print(f"    -> 没有找到匹配的行")
        
    except Exception as e:
        print(f"    错误: {str(e)}")

def extract_meth_data(features, input_file, output_file, meth_type):
    """从甲基化bedGraph文件中提取特征数据(基于行号,无表头)"""
    print(f"\n  处理{meth_type}甲基化文件: {os.path.basename(input_file)}")
    
    # 筛选出对应类型的甲基化特征
    meth_features = {f for f in features if f.startswith(meth_type + '_')}
    if not meth_features:
        print(f"    没有找到{meth_type}特征")
        return
    
    print(f"    需要查找 {len(meth_features)} 个{meth_type}特征")
    
    try:
        # 读取bedGraph文件（制表符分隔，无表头）
        # 先读取所有列
        df = pd.read_csv(input_file, sep='\t', header=None)
        
        # 选择需要的列：chr(第1列), start(第2列), end(第3列), CV(最后一列)
        # 列索引从0开始
        df_selected = pd.DataFrame({
            'chr': df.iloc[:, 0],      # 第1列：染色体
            'start': df.iloc[:, 1],    # 第2列：起始位置
            'end': df.iloc[:, 2],      # 第3列：结束位置
            'CV': df.iloc[:, -1]       # 最后一列：CV值
        })
        
        # 提取行号（特征数字直接对应行号，从1开始）
        rows_to_extract = []
        feature_names = []
        
        for feature in sorted(meth_features):
            _, pos = parse_feature_name(feature)
            if pos is not None:
                # 甲基化文件无表头，行号从1开始，pandas索引从0开始
                # 所以实际索引 = 行号 - 1
                row_index = pos - 1
                if 0 <= row_index < len(df_selected):
                    rows_to_extract.append(row_index)
                    feature_names.append(feature)
        
        # 提取对应行
        if rows_to_extract:
            matched_df = df_selected.iloc[rows_to_extract].copy()
            matched_df.insert(0, 'Feature', feature_names)
            
            # 统一chr列为Chr2
            matched_df['chr'] = 'Chr2'
            
            # 按照start列（第3列）升序排序
            matched_df = matched_df.sort_values('start')
            
            # 保存结果（CSV格式，逗号分隔）
            matched_df.to_csv(output_file, index=False)
            print(f"    -> 找到 {len(matched_df)} 行数据")
            print(f"    -> 保存到: {output_file}")
        else:
            print(f"    -> 没有找到匹配的行")
        
    except Exception as e:
        print(f"    错误: {str(e)}")

def process_feature_category(category_name, feature_file, data_base_dir, output_dir):
    """处理单个特征类别"""
    print(f"\n{'='*60}")
    print(f"处理特征类别: {category_name}")
    print(f"{'='*60}")
    
    # 读取特征文件
    features = read_feature_csv(feature_file)
    if not features:
        print(f"  跳过 {category_name}（无特征）")
        return
    
    # 创建输出子目录
    category_output_dir = os.path.join(output_dir, category_name)
    os.makedirs(category_output_dir, exist_ok=True)
    
    # 1. 处理SNP数据
    snp_input = os.path.join(data_base_dir, "snp/filtered_snp_frequencies.csv")
    snp_output = os.path.join(category_output_dir, f"{category_name}_snp_data.csv")
    if os.path.exists(snp_input):
        extract_snp_data(features, snp_input, snp_output)
    else:
        print(f"\n  警告: SNP文件不存在 - {snp_input}")
    
    # 2. 处理CG甲基化数据
    cg_input = os.path.join(data_base_dir, "meth/CG.CV.filter.bedGraph")
    cg_output = os.path.join(category_output_dir, f"{category_name}_CG_data.csv")
    if os.path.exists(cg_input):
        extract_meth_data(features, cg_input, cg_output, "CG")
    else:
        print(f"\n  警告: CG文件不存在 - {cg_input}")
    
    # 3. 处理CHG甲基化数据
    chg_input = os.path.join(data_base_dir, "meth/CHG.CV.filter.bedGraph")
    chg_output = os.path.join(category_output_dir, f"{category_name}_CHG_data.csv")
    if os.path.exists(chg_input):
        extract_meth_data(features, chg_input, chg_output, "CHG")
    else:
        print(f"\n  警告: CHG文件不存在 - {chg_input}")
    
    # 4. 处理CHH甲基化数据
    chh_input = os.path.join(data_base_dir, "meth/CHH.CV.filter.bedGraph")
    chh_output = os.path.join(category_output_dir, f"{category_name}_CHH_data.csv")
    if os.path.exists(chh_input):
        extract_meth_data(features, chh_input, chh_output, "CHH")
    else:
        print(f"\n  警告: CHH文件不存在 - {chh_input}")

def main():
    # 基础路径
    feature_dir = "../../list/xgboot/TPM_4.5/feature_analysis"
    data_base_dir = "../../list"
    output_dir = "../../list/xgboot/TPM_4.5/feature_data_extraction"
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("特征数据提取")
    print("=" * 60)
    
    # 定义需要处理的特征文件
    feature_categories = {
        "Sexpr_only": os.path.join(feature_dir, "Sexpr_only.csv"),
        "Ssplice_only": os.path.join(feature_dir, "Ssplice_only.csv"),
        "Scoupled_specific": os.path.join(feature_dir, "Scoupled_specific.csv")
    }
    
    # 处理每个特征类别
    for category_name, feature_file in feature_categories.items():
        if os.path.exists(feature_file):
            process_feature_category(category_name, feature_file, 
                                    data_base_dir, output_dir)
        else:
            print(f"\n警告: 特征文件不存在 - {feature_file}")
    
    # 生成总结报告
    print(f"\n{'='*60}")
    print("生成总结报告")
    print(f"{'='*60}")
    
    report = f"""
特征数据提取报告
{'=' * 60}

输出目录: {output_dir}

提取的特征类别:
1. Sexpr_only (表达量主导)
   - {output_dir}/Sexpr_only/Sexpr_only_snp_data.csv
   - {output_dir}/Sexpr_only/Sexpr_only_CG_data.csv
   - {output_dir}/Sexpr_only/Sexpr_only_CHG_data.csv
   - {output_dir}/Sexpr_only/Sexpr_only_CHH_data.csv

2. Ssplice_only (剪接主导)
   - {output_dir}/Ssplice_only/Ssplice_only_snp_data.csv
   - {output_dir}/Ssplice_only/Ssplice_only_CG_data.csv
   - {output_dir}/Ssplice_only/Ssplice_only_CHG_data.csv
   - {output_dir}/Ssplice_only/Ssplice_only_CHH_data.csv

3. Scoupled_specific (耦合)
   - {output_dir}/Scoupled_specific/Scoupled_specific_snp_data.csv
   - {output_dir}/Scoupled_specific/Scoupled_specific_CG_data.csv
   - {output_dir}/Scoupled_specific/Scoupled_specific_CHG_data.csv
   - {output_dir}/Scoupled_specific/Scoupled_specific_CHH_data.csv

说明:
- SNP文件有表头，特征snp_N对应第N+1行（包含表头）
- 甲基化文件无表头，特征CG_N/CHG_N/CHH_N对应第N行
- 所有文件行号从1开始计数

{'=' * 60}
"""
    
    # 保存报告
    report_file = os.path.join(output_dir, "extraction_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(report)
    print(f"\n报告已保存到: {report_file}")
    print("\n数据提取完成！")

if __name__ == "__main__":
    main()
