from sklearn.datasets import fetch_california_housing
from sklearn.model_selection import train_test_split
from sklearn import metrics
import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost import plot_importance
import matplotlib.pyplot as plt
import scipy.stats as stats

# 加载加州房价数据集
housing = fetch_california_housing()
X, y = housing.data, housing.target

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# 定义XGBoost回归模型
xgb_regressor = xgb.XGBRegressor(
    objective='reg:squarederror',  # 回归问题的目标函数
    max_depth=5,                  # 树的最大深度
    learning_rate=0.1,            # 学习率
    n_estimators=100,             # 树的个数
    subsample=0.8,                # 训练每棵树时使用的样本比例
    colsample_bytree=0.8          # 训练每棵树时使用的特征比例
)

# 训练模型
xgb_regressor.fit(X_train, y_train)

# 预测测试集
y_pred = xgb_regressor.predict(X_test)

# 计算均方误差
mse = metrics.mean_squared_error(y_test, y_pred)
#print(f"Mean Squared Error: {mse}")

# 也可以计算均方根误差（Root Mean Squared Error, RMSE）
rmse = np.sqrt(mse)
#print(f"Root Mean Squared Error: {rmse}")

# 计算平均绝对误差（Mean Absolute Error, MAE）
mae = metrics.mean_absolute_error(y_test, y_pred)
#print(f"Mean Absolute Error: {mae}")

# 计算决定系数（R²）
r2 = metrics.r2_score(y_test, y_pred)
#print(f"R²: {r2}")

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
plt.savefig('feature_importance.png')  # 保存特征重要性图
plt.show()


# 计算皮尔逊相关系数
pearson_corr, _ = stats.pearsonr(y_test, y_pred)
print(f"Pearson Correlation Coefficient: {pearson_corr}")

# 交叉验证
####交叉验证####
params = {
    'max_depth': 10,
    'eta': 0.05,
    'objective': 'reg:squarederror',
    'eval_metric': 'rmse'  # 注意：cv返回的是mse，但这里指定rmse是为了让xgboost内部监控
}
dtrain = xgb.DMatrix(X_train, label=y_train)
cv_results = xgb.cv(params, dtrain, num_boost_round=200, nfold=10, metrics={'rmse'}, early_stopping_rounds=10, seed=42)
best_nrounds = cv_results.shape[0] - 1


# 用最佳迭代次数训练模型
bst = xgb.train(params, dtrain, num_boost_round=best_nrounds)
dtest = xgb.DMatrix(X_test)
y_pred_cv = bst.predict(dtest)


# 计算均方误差
mse = metrics.mean_squared_error(y_test, y_pred_cv)
#print(f"Mean Squared Error: {mse}")

# 也可以计算均方根误差（Root Mean Squared Error, RMSE）
rmse = np.sqrt(mse)
#print(f"Root Mean Squared Error: {rmse}")

# 计算平均绝对误差（Mean Absolute Error, MAE）
mae = metrics.mean_absolute_error(y_test, y_pred_cv)
#print(f"Mean Absolute Error: {mae}")

# 计算决定系数（R²）
r2 = metrics.r2_score(y_test, y_pred_cv)
#print(f"R²: {r2}")

print("交叉验证后均方误差 (MSE):", mse)
print("交叉验证后均方根误差 (RMSE):", rmse)
print("交叉验证后平均绝对误差 (MAE):", mae)
print("交叉验证后拟合优度 (R-squared):", r2)

# 可视化特征重要性
plt.figure()
plot_importance(xgb_regressor, ax=plt.gca())
plt.title('Feature Importance')
plt.savefig('feature_importance_cv.png')  # 保存特征重要性图
plt.show()

# 计算皮尔逊相关系数
pearson_corr_cv, _ = stats.pearsonr(y_test, y_pred_cv)
print(f"Pearson Correlation Coefficient after CV: {pearson_corr_cv}")

# 绘制散点图
plt.figure()
plt.scatter(y_test, y_pred)

# 添加一些图表装饰
plt.xlabel('True Values (y_test)')
plt.ylabel('Predicted Values (y_pred)')
plt.title('Scatter Plot of True vs Predicted Values')
plt.grid(True)

# 可选：添加一条理想情况下的对角线（即完全准确的预测）
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=4)

plt.savefig('scatter_plot.png')  # 保存散点图
# 显示图表
plt.show()

# 保存模型
bst.save_model('my_model.json')
