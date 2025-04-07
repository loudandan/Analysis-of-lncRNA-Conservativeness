import os
import numpy as np
import pandas as pd
# from joblib import load

# # 加载训练好的模型
# model = load('random_forest_model.joblib')  # 假设模型已保存


import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix

# 加载数据
features = np.load('combined_features.npy')  # 形状 (79611, 12)
labels = np.load('combined_labels.npy').ravel()  # 转换为一维数组 (79611,)

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(
    features, labels, 
    test_size=0.2, 
    random_state=42, 
    stratify=labels
)

# 初始化随机森林模型
rf = RandomForestClassifier(
    n_estimators=100,
    random_state=42,
    n_jobs=-1  # 使用所有CPU核心加速训练
)

# 训练模型
rf.fit(X_train, y_train)

# 预测测试集
y_pred = rf.predict(X_test)

# 评估模型
print("测试集准确率:", accuracy_score(y_test, y_pred))
print("\n分类报告:")
print(classification_report(y_test, y_pred))
print("\n混淆矩阵:")
print(confusion_matrix(y_test, y_pred))

# 特征重要性
print("\n特征重要性:")
for i, importance in enumerate(rf.feature_importances_):
    print(f"特征 {i}: {importance:.4f}")

# 定义种类列表和路径
kinds = ["anwen", "bm", "cngd", "ls", "mifeng", "sf", "yachong"]
root_dir = './预测集/prediction/'
output_dir = './prediction_results'
os.makedirs(output_dir, exist_ok=True)

# 定义特征列索引（根据文件格式调整）
feature_columns = [
    'mu1', 'md1', 'mf1', 'proportionscoreu1', 'proportionscored1', 'proportionscoref1',
    'mu2', 'md2', 'mf2', 'proportionscoreu2', 'proportionscored2', 'proportionscoref2'
]

for kind in kinds:
    # 构建文件路径
    file_path = os.path.join(root_dir, f'nl_{kind}_predict_pred_nl_insect_feature.csv')
    print(f"处理文件: {file_path}")
    try:
        # 读取数据文件
        df = pd.read_csv(file_path, sep=',', header=0)
        
        # 提取特征数据
        X_pred = df[feature_columns].to_numpy()
        
        # 检查数据维度
        if X_pred.shape[1] != 12:
            raise ValueError(f"文件 {kind} 特征维度错误，预期12列，实际得到{X_pred.shape[1]}列")
        
        # 执行预测
        predictions = rf.predict(X_pred)
        
        # 构建结果DataFrame
        result_df = pd.DataFrame({
            'nl_gene1': df['nl_gene1'],
            'insect_gene1': df['insect_gene1'],
            'class': predictions
        })
        
        # 保存结果
        output_path = os.path.join(output_dir, f'{kind}_predictions.tsv')
        result_df.to_csv(output_path, sep='\t', index=False)
        print(f"成功处理 {kind}，结果已保存至 {output_path}")
        
    except Exception as e:
        print(f"处理 {kind} 时发生错误: {str(e)}")
        continue

print("\n所有文件处理完成！")