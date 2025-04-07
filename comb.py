# # 读取对应关系文件
# mappings = []
# with open('nl_anwen_final.tsv') as f:
#     next(f)  # 跳过标题行
#     for line in f:
#         parts = line.strip().split('\t')
#         if len(parts) >= 3:
#             nl_id = parts[1]  # 第二列
#             ag_id = parts[2]  # 第三列
#             # print(nl_id,ag_id)
#             mappings.append((nl_id, ag_id))




# # 读取test1.bed数据
# test1_data = {}
# with open('nlgenome.bed') as f:
#     for line in f:
#         cols = line.strip().split('\t')
#         if len(cols) >= 5:
#             key = cols[4]
#             test1_data[key] = cols[:4]  # 存储前四列
#             print(cols)

# # 读取test2.bed数据
# test2_data = {}
# with open('anwengenome.bed') as f:
#     for line in f:
#         cols = line.strip().split('\t')
#         if len(cols) >= 5:
#             # print(cols)
#             key = cols[4]
#             test2_data[key] = cols[:4]  # 存储前四列

# # 生成合并结果
# result = []
# for nl_id, ag_id in mappings:
#     if nl_id in test1_data and ag_id in test2_data:
#         combined = test1_data[nl_id] + test2_data[ag_id]
#         result.append('\t'.join(combined))

# # # # 打印结果（或写入文件）
# # print('\n'.join(result))

# with open('comb_res.txt','w') as f:
#     f.write('\n'.join(result))






import argparse

def read_mappings(mapping_file):
    """读取ID映射关系文件"""
    mappings = []
    with open(mapping_file) as f:
        next(f)  # 跳过标题行
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                mappings.append((parts[1], parts[2]))  # (nl_id, ag_id)
                # print(parts[1], parts[2])

    return mappings

def read_bed_file(bed_file):
    """读取BED格式文件，返回{ID: 前四列数据}的字典"""
    data = {}
    with open(bed_file) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 5:
                data[cols[4]] = cols[:4]  # 使用第5列作为键
                # print(cols)
    return data

def main(mapping_file, bed1_file, bed2_file, output_file):
    """主处理函数"""
    # 读取数据
    mappings = read_mappings(mapping_file)
    bed1_data = read_bed_file(bed1_file)
    bed2_data = read_bed_file(bed2_file)
    
    # 合并匹配结果
    result = []
    for nl_id, ag_id in mappings:
        if nl_id in bed1_data and ag_id in bed2_data:
            combined = bed1_data[nl_id] + bed2_data[ag_id]
            result.append('\t'.join(combined))
            print(combined)
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        f.write('\n'.join(result))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='合并BED文件数据工具')
    parser.add_argument('mapping', help='映射关系文件(TSV格式)')
    parser.add_argument('bed1', help='第一个BED文件')
    parser.add_argument('bed2', help='第二个BED文件')
    parser.add_argument('output', help='输出文件路径')
    
    args = parser.parse_args()
    
    main(
        mapping_file=args.mapping,
        bed1_file=args.bed1,
        bed2_file=args.bed2,
        output_file=args.output
    )