#!/usr/bin/env python3
import csv
import sys

def process_csv(input_file, output_file):
    """
    处理CSV文件：
    1. 提取第1列(name)和第2列(MCC)
    2. 去除引号
    3. 输出为制表符分隔
    """
    data_rows = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # 跳过标题行
        
        for row in reader:
            if len(row) >= 9:  # 确保有足够的列
                name = row[0].strip('"')  # 第1列（索引0），去除引号
                mcc = row[1].strip('"')  # 第2列（索引1），去除引号
                
                data_rows.append([name, mcc])
    
    # 按照第2列（MCC值）从高到低排序
    data_rows.sort(key=lambda x: float(x[1]), reverse=True)
    
    # 写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入数据行（不写标题行）
        for row in data_rows:
            f.write(f"{row[0]}\t{row[1]}\n")
    
    print(f"处理完成！输出文件：{output_file}")
    print(f"共处理 {len(data_rows)} 行数据")

# total
if __name__ == "__main__":
    input_file = "total/total.MCC.csv"
    output_file = "total/total.gene_id.txt"
    
    try:
        process_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

if __name__ == "__main__":
    input_file = "total/total.negative.MCC.csv"
    output_file = "total/total.negative.gene_id.txt"
    
    try:
        process_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

# 10C
    input_file_10C = "10C/10C.MCC.csv"
    output_file_10C = "10C/10C.gene_id.txt"
    
    try:
        process_csv(input_file_10C, output_file_10C)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file_10C}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

    input_file_10C = "10C/10C.negative.MCC.csv"
    output_file_10C = "10C/10C.negative.gene_id.txt"
    
    try:
        process_csv(input_file_10C, output_file_10C)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file_10C}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

# 16C
    input_file_16C = "16C/16C.MCC.csv"
    output_file_16C = "16C/16C.gene_id.txt"
    
    try:
        process_csv(input_file_16C, output_file_16C)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file_16C}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

    input_file_16C = "16C/16C.negative.MCC.csv"
    output_file_16C = "16C/16C.negative.gene_id.txt"
    
    try:
        process_csv(input_file_16C, output_file_16C)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file_16C}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

# 22C
    input_file_22C = "22C/22C.MCC.csv"
    output_file_22C = "22C/22C.gene_id.txt"

    try:
        process_csv(input_file_22C, output_file_22C)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file_22C}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

    input_file_22C = "22C/22C.negative.MCC.csv"
    output_file_22C = "22C/22C.negative.gene_id.txt"

    try:
        process_csv(input_file_22C, output_file_22C)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file_22C}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")
