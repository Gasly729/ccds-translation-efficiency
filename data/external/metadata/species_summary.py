#!/usr/bin/env python3
"""
统计 metadata.csv 中的物种分布
"""

import csv
from collections import Counter

def main():
    metadata_file = "/home/xrx/my_projects/joker/26.3/download/metadata/metadata.csv"
    output_file = "/home/xrx/my_projects/joker/26.3/download/metadata/species_summary.csv"
    
    species_counts = Counter()
    total_records = 0
    
    with open(metadata_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        
        # 跳过第一行（注释行）
        next(reader)
        
        # 读取表头
        headers = next(reader)
        organism_col = headers.index("organism")
        
        for row_num, row in enumerate(reader, start=3):  # 从第3行开始计数
            if len(row) <= organism_col:
                continue
                
            organism = row[organism_col].strip()
            if organism:
                species_counts[organism] += 1
                total_records += 1
    
    # 写入统计结果
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['organism', 'count', 'percentage'])
        
        for organism, count in species_counts.most_common():
            percentage = (count / total_records) * 100
            writer.writerow([organism, count, f"{percentage:.2f}%"])
    
    print(f"物种统计完成！")
    print(f"总记录数: {total_records}")
    print(f"物种数: {len(species_counts)}")
    print(f"输出文件: {output_file}")
    
    # 打印前10个物种
    print("\n物种分布前10:")
    for i, (organism, count) in enumerate(species_counts.most_common(10), 1):
        percentage = (count / total_records) * 100
        print(f"{i:2d}. {organism}: {count} ({percentage:.2f}%)")

if __name__ == "__main__":
    main()
