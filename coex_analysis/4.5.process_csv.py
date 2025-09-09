#!/usr/bin/env python3
import csv
import sys

# all
def process_csv(input_file, output_file):
    """
    处理CSV文件：
    1. 提取第9列(name)和第6列(Degree)
    2. 去除引号
    3. 输出为制表符分隔
    """
    data_rows = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # 跳过标题行
        
        for row in reader:
            if len(row) >= 9:  # 确保有足够的列
                name = row[8].strip('"')  # 第9列（索引8），去除引号
                degree = row[5].strip('"')  # 第6列（索引5），去除引号
                
                data_rows.append([name, degree])
    
    # 写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入数据行（不写标题行）
        for row in data_rows:
            f.write(f"{row[0]}\t{row[1]}\n")
    
    print(f"处理完成！输出文件：{output_file}")
    print(f"共处理 {len(data_rows)} 行数据")

if __name__ == "__main__":
    input_file = "total/coex-brown-salmon.node.csv"
    output_file = "total/all.gene_id.txt"
    
    try:
        process_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

# A
def process_csv(input_file, output_file):
    """
    处理CSV文件：
    1. 提取第9列(name)和第6列(Degree)
    2. 去除引号
    3. 输出为制表符分隔
    """
    data_rows = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # 跳过标题行
        
        for row in reader:
            if len(row) >= 9:  # 确保有足够的列
                name = row[8].strip('"')  # 第9列（索引8），去除引号
                degree = row[5].strip('"')  # 第6列（索引5），去除引号
                
                data_rows.append([name, degree])
    
    # 写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入数据行（不写标题行）
        for row in data_rows:
            f.write(f"{row[0]}\t{row[1]}\n")
    
    print(f"处理完成！输出文件：{output_file}")
    print(f"共处理 {len(data_rows)} 行数据")

if __name__ == "__main__":
    input_file = "a/a.coex-yellow.node.csv"
    output_file = "a/a.gene_id.txt"
    
    try:
        process_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")

# B
def process_csv(input_file, output_file):
    """
    处理CSV文件：
    1. 提取第9列(name)和第6列(Degree)
    2. 去除引号
    3. 输出为制表符分隔
    """
    data_rows = []
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)  # 跳过标题行
        
        for row in reader:
            if len(row) >= 9:  # 确保有足够的列
                name = row[8].strip('"')  # 第9列（索引8），去除引号
                degree = row[5].strip('"')  # 第6列（索引5），去除引号
                
                data_rows.append([name, degree])
    
    # 写入输出文件
    with open(output_file, 'w', encoding='utf-8') as f:
        # 写入数据行（不写标题行）
        for row in data_rows:
            f.write(f"{row[0]}\t{row[1]}\n")
    
    print(f"处理完成！输出文件：{output_file}")
    print(f"共处理 {len(data_rows)} 行数据")

if __name__ == "__main__":
    input_file = "b/b.coex-blue-pink.node.csv"
    output_file = "b/b.gene_id.txt"

    try:
        process_csv(input_file, output_file)
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {input_file}")
    except Exception as e:
        print(f"处理过程中出现错误：{e}")
