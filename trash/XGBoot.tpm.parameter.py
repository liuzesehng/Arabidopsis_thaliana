import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import metrics
import xgboost as xgb
from xgboost import plot_importance
import matplotlib.pyplot as plt
import scipy.stats as stats

# 数据读取
df = pd.read_csv('../RCA/RCA.tsv', sep='\t')
df.head()

# 缺失值处理
df.isnull().sum()

# 划分训练集、验证集和测试集
X_temp, X_test, y_temp, y_test = train_test_split(df.iloc[:, [0, 1, 2]], df['TPM'], test_size=0.2, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.125, random_state=42)  # 0.125 x 0.8 = 0.1

# 输出数据集的大小
print(f"训练集维度: {X_train.shape}")
print(f"验证集维度: {X_val.shape}")
print(f"测试集维度: {X_test.shape}")

# 参数范围设置
param_grid = {
    'max_depth': [3, 4, 5, 6, 7, 8, 9, 10],
    'learning_rate': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
    'n_estimators': [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],
    'subsample': [0.7, 0.8, 0.9],
    'colsample_bytree': [0.7,0.8, 0.9]
}

# 循环遍历参数组合并训练模型
best_rmse = float("inf")
best_params = None

for max_depth in param_grid['max_depth']:
    for learning_rate in param_grid['learning_rate']:
        for n_estimators in param_grid['n_estimators']:
            for subsample in param_grid['subsample']:
                for colsample_bytree in param_grid['colsample_bytree']:
                    params = {
                        'objective': 'reg:squarederror',
                        'max_depth': max_depth,
                        'learning_rate': learning_rate,
                        'n_estimators': n_estimators,
                        'subsample': subsample,
                        'colsample_bytree': colsample_bytree,
                        'early_stopping_rounds': 10,
                        'eval_metric': 'rmse'
                    }
                    
                    xgb_regressor = xgb.XGBRegressor(**params)
                    
                    # 训练模型
                    xgb_regressor.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)
                    
                    # 预测测试集
                    y_test_pred = xgb_regressor.predict(X_test)
                    
                    # 计算RMSE
                    rmse = np.sqrt(metrics.mean_squared_error(y_test, y_test_pred))
                    
                    # 打印当前参数组合和RMSE
                    print(f"Params: {params}, RMSE: {rmse}")
                    
                    # 记录最优参数组合
                    if rmse < best_rmse:
                        best_rmse = rmse
                        best_params = params

print(f"Best params: {best_params}, Best RMSE: {best_rmse}")

# 用最佳参数组合重新训练模型
xgb_regressor = xgb.XGBRegressor(**best_params)
xgb_regressor.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)

# 预测测试集
y_test_pred = xgb_regressor.predict(X_test)

# 计算测试集评估指标
mse = metrics.mean_squared_error(y_test, y_test_pred)
rmse = np.sqrt(mse)
mae = metrics.mean_absolute_error(y_test, y_test_pred)
r2 = metrics.r2_score(y_test, y_test_pred)

print("测试集评估指标")
print("均方误差 (MSE):", mse)
print("均方根误差 (RMSE):", rmse)
print("平均绝对误差 (MAE):", mae)
print("拟合优度 (R-squared):", r2)

# 可视化特征重要性
plt.figure()
plot_importance(xgb_regressor, ax=plt.gca())
plt.title('Feature Importance')
plt.xlabel('F-Score')  # 设置横坐标轴标题
plt.ylabel('Features')  # 设置纵坐标轴标题
plt.savefig('feature_importance_TPM_1.0.png')  # 保存特征重要性图
plt.show()

# 计算皮尔逊相关系数
pearson_corr, _ = stats.pearsonr(y_test, y_test_pred)
print(f"Pearson Correlation Coefficient: {pearson_corr}")

# 绘制散点图
plt.figure()
plt.scatter(y_test, y_test_pred)

# 添加一些图表装饰
plt.xlabel('True Values (y_test)')
plt.ylabel('Predicted Values (y_test_pred)')
plt.title('Scatter Plot of True vs Predicted Values')
plt.grid(True)

# 可选：添加一条理想情况下的对角线（即完全准确的预测）
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=4)

plt.savefig('scatter_plot_TPM_1.0.png')  # 保存散点图
# 显示图表
plt.show()

# 保存模型
xgb_regressor.save_model('my_model_TPM_1.0.json')
