#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
特征集合分析脚本
从XGBoost优化结果中提取特征，进行集合运算并绘制韦恩图
"""

import re
import pandas as pd
import os

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

def save_features_to_csv(feature_set, output_file, column_name='Feature'):
    """将特征集合保存为CSV文件"""
    df = pd.DataFrame(sorted(list(feature_set)), columns=[column_name])
    df.to_csv(output_file, index=False)
    print(f"  保存CSV: {output_file} ({len(feature_set)} 个特征)")

def main():
    # 基础路径
    base_path = "../../list/xgboot/TPM_4.5"
    output_dir = "../../list/xgboot/TPM_4.5/feature_analysis"
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("特征集合分析")
    print("=" * 60)
    
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
    
    # Sexpr-core = SA ∩ SB ∩ ST (表达量核心)
    Sexpr_core = SA & SB & ST
    print(f"  Sexpr-core (SA ∩ SB ∩ ST): {len(Sexpr_core)} 个特征")
    
    # Ssplice-core = SPα ∩ SPβ (剪接核心)
    Ssplice_core = SPa & SPb
    print(f"  Ssplice-core (SPα ∩ SPβ): {len(Ssplice_core)} 个特征")
    
    # ========================================
    # 3. 计算特异性集合
    # ========================================
    print("\n[步骤 3] 计算特异性集合...")
    
    # Sexpr-only = Sexpr-core \ Ssplice-core (表达量主导)
    Sexpr_only = Sexpr_core - Ssplice_core
    print(f"  Sexpr-only (表达量主导): {len(Sexpr_only)} 个特征")
    
    # Ssplice-only = Ssplice-core \ Sexpr-core (剪接主导)
    Ssplice_only = Ssplice_core - Sexpr_core
    print(f"  Ssplice-only (剪接主导): {len(Ssplice_only)} 个特征")
    
    # Scoupled-specific = Sexpr-core ∩ Ssplice-core (耦合)
    Scoupled_specific = Sexpr_core & Ssplice_core
    print(f"  Scoupled-specific (耦合): {len(Scoupled_specific)} 个特征")
    
    # ========================================
    # 4. 保存CSV文件
    # ========================================
    print("\n[步骤 4] 保存CSV文件...")
    
    # 保存核心集合
    save_features_to_csv(Sexpr_core, 
                        f"{output_dir}/Sexpr_core.csv",
                        "Expression_Core_Feature")
    
    save_features_to_csv(Ssplice_core,
                        f"{output_dir}/Ssplice_core.csv",
                        "Splicing_Core_Feature")
    
    # 保存特异性集合
    save_features_to_csv(Sexpr_only,
                        f"{output_dir}/Sexpr_only.csv",
                        "Expression_Only_Feature")
    
    save_features_to_csv(Ssplice_only,
                        f"{output_dir}/Ssplice_only.csv",
                        "Splicing_Only_Feature")
    
    save_features_to_csv(Scoupled_specific,
                        f"{output_dir}/Scoupled_specific.csv",
                        "Coupled_Feature")
    
    # ========================================
    # 5. 生成统计报告
    # ========================================
    print("\n[步骤 5] 生成统计报告...")
    
    report = f"""
特征集合分析报告
{'=' * 60}

1. 原始特征集合:
   - SA (A_TPM): {len(SA)} 个特征
   - SB (B_TPM): {len(SB)} 个特征
   - ST (total_TPM): {len(ST)} 个特征
   - SPα (A_per_TPM): {len(SPa)} 个特征
   - SPβ (B_per_TPM): {len(SPb)} 个特征

2. 核心集合:
   - Sexpr-core (SA ∩ SB ∩ ST): {len(Sexpr_core)} 个特征
   - Ssplice-core (SPα ∩ SPβ): {len(Ssplice_core)} 个特征

3. 特异性集合:
   - Sexpr-only (表达量主导): {len(Sexpr_only)} 个特征
     = Sexpr-core \\ Ssplice-core
   
   - Ssplice-only (剪接主导): {len(Ssplice_only)} 个特征
     = Ssplice-core \\ Sexpr-core
   
   - Scoupled-specific (耦合): {len(Scoupled_specific)} 个特征
     = Sexpr-core ∩ Ssplice-core

4. 验证:
   - 表达核心 - 耦合 = 表达独有: {len(Sexpr_core)} - {len(Scoupled_specific)} = {len(Sexpr_core) - len(Scoupled_specific)}
   - 实际表达独有: {len(Sexpr_only)}
   - 一致性检查: {'通过' if len(Sexpr_only) == len(Sexpr_core) - len(Scoupled_specific) else '失败'}

5. 输出文件:
   CSV文件:
   - {output_dir}/Sexpr_core.csv
   - {output_dir}/Ssplice_core.csv
   - {output_dir}/Sexpr_only.csv
   - {output_dir}/Ssplice_only.csv
   - {output_dir}/Scoupled_specific.csv

{'=' * 60}
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
