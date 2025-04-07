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