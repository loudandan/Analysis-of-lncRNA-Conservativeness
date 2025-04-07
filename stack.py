import os
import numpy as np
import pandas as pd

# 定义种类列表
kinds = ["anwen", "bm", "cngd", "ls", "mifeng", "sf", "yachong"]

# 初始化数据存储列表
data_list = []
label_list = []

# 定义需要提取的列名
required_columns = [
    'mu1', 'md1', 'mf1', 'proportionscoreu1', 'proportionscored1', 'proportionscoref1',
    'mu2', 'md2', 'mf2', 'proportionscoreu2', 'proportionscored2', 'proportionscoref2'
]

col_ranges = [
        slice(2, 8),   # 第3-8列（索引2到7）
        slice(10, 16)  # 第11-16列（索引10到15）
    ]


for kind in kinds:
    # 构建文件路径
    root_dir = os.path.join('./训练集', f'nl_{kind}')
    pos_path = os.path.join(root_dir, f'{kind}train_nl_insect_pos_feature.csv')
    neg_path = os.path.join(root_dir, f'{kind}train_nl_insect_neg_feature.csv')
    
    # 处理正样本文件
    if os.path.exists(pos_path):
        try:
            df = pd.read_csv(pos_path, sep=',')
            # print(df.iloc[:, 0] )
            data = df[required_columns].values.astype(np.float32)

            # # 合并目标列
            # selected_cols = np.hstack([
            #     df.iloc[:, col_range].values 
            #     for col_range in col_ranges
            # ])
            # # print([
            # #     df.iloc[:, col_range].values 
            # #     for col_range in col_ranges
            # # ])
            # data = selected_cols.astype(np.float32)
            print(data)
            data_list.append(data)
            label_list.append(np.ones((len(data), 1), dtype=np.int8))
        except Exception as e:
            print(f"Error processing {pos_path}: {str(e)}")
    
    # 处理负样本文件
    if os.path.exists(neg_path):
        try:
            df = pd.read_csv(neg_path, sep=',')
            data = df[required_columns].values.astype(np.float32)
            data_list.append(data)
            label_list.append(np.zeros((len(data), 1), dtype=np.int8))
        except Exception as e:
            print(f"Error processing {neg_path}: {str(e)}")

# 合并所有数据
final_data = np.vstack(data_list) if data_list else np.array([])
final_labels = np.vstack(label_list) if label_list else np.array([])

# 保存结果
np.save('combined_features.npy', final_data)
np.save('combined_labels.npy', final_labels)

print(f"数据合并完成，最终维度：{final_data.shape}")
print(f"标签维度：{final_labels.shape}")