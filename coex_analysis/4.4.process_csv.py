#!/usr/bin/env python3
import csv
import sys

def process_csv(input_file, output_file):
    """
    处理CSV文件：
    1. 提取第2列(name)和最后一列(weight)
    2. 去除引号
    3. 将name列按" (interacts with) "分割
    4. AT2G39730放第一列，其他放第二列
    5. 按weight降序排序
    6. 输出为制表符分隔
    """
    data_rows = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # 跳过标题行
        
        for row in reader:
            if len(row) >= 6:  # 确保有足够的列
                name = row[1].strip('"')  # 第2列，去除引号
                weight = float(row[5].strip('"'))  # 最后一列，去除引号并转换为数字
                
                # 按" (interacts with) "分割name列
                if " (interacts with) " in name:
                    parts = name.split(" (interacts with) ")
                    if len(parts) == 2:
                        gene1, gene2 = parts[0], parts[1]
                        
                        # AT2G39730放在第一列
                        if gene1 == "AT2G39730":
                            first_col = gene1
                            second_col = gene2
                        elif gene2 == "AT2G39730":
                            first_col = "AT2G39730"
                            second_col = gene1
                        else:
                            # 如果都不是AT2G39730，保持原顺序
                            first_col = gene1
                            second_col = gene2
                        
                        data_rows.append([first_col, second_col, weight])
    
    # 按weight降序排序
    data_rows.sort(key=lambda x: x[2], reverse=True)
    
    # 写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入标题行
        f.write("name\tshared_name\tweight\n")
        
        # 写入数据行
        for row in data_rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
    
    print(f"处理完成！输出文件：{output_file}")
    print(f"共处理 {len(data_rows)} 行数据")

if __name__ == "__main__":
    input_file = "total/CytoscapeInput-edges-rca.csv"
    output_file = "total/gene_id.tsv"
    
    try:
        process_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")
