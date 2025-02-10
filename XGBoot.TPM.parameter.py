import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
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

# 划分训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(df.iloc[:,[0, 1, 2]], df['TPM'], test_size=0.2, random_state=42)

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
                    }
                    
                    xgb_regressor = xgb.XGBRegressor(**params)
                    
                    # 训练模型
                    xgb_regressor.fit(X_train, y_train)
                    
                    # 预测测试集
                    y_pred = xgb_regressor.predict(X_test)
                    
                    # 计算RMSE
                    rmse = np.sqrt(metrics.mean_squared_error(y_test, y_pred))
                    
                    # 打印当前参数组合和RMSE
                    print(f"Params: {params}, RMSE: {rmse}")
                    
                    # 记录最优参数组合
                    if rmse < best_rmse:
                        best_rmse = rmse
                        best_params = params

print(f"Best params: {best_params}, Best RMSE: {best_rmse}")

# 用最佳参数组合重新训练模型
xgb_regressor = xgb.XGBRegressor(**best_params)
xgb_regressor.fit(X_train, y_train)


# 预测测试集
y_pred = xgb_regressor.predict(X_test)

# 计算评估指标
mse = metrics.mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
mae = metrics.mean_absolute_error(y_test, y_pred)
r2 = metrics.r2_score(y_test, y_pred)

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
pearson_corr, _ = stats.pearsonr(y_test, y_pred)
print(f"Pearson Correlation Coefficient: {pearson_corr}")


# 交叉验证：
param_grid = {
    'max_depth': [3, 4, 5, 6, 7, 8, 9, 10],
    'learning_rate': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
    'subsample': [0.7, 0.8, 0.9],
    'colsample_bytree': [0.7, 0.8, 0.9]
}


# 循环遍历参数组合并训练模型
best_rmse = float("inf")
best_params = None

for max_depth in param_grid['max_depth']:
    for learning_rate in param_grid['learning_rate']:
        params = {
            'objective': 'reg:squarederror',
            'max_depth': max_depth,
            'learning_rate': learning_rate,
            'eval_metric': 'rmse'  # 注意：cv返回的是mse，但这里指定rmse是为了让xgboost内部监控
        }
              
        # 转换数据为DMatrix格式
        dtrain = xgb.DMatrix(X_train, label=y_train)
        cv_results = xgb.cv(params, dtrain, num_boost_round=1000, nfold=10, metrics={'rmse'}, early_stopping_rounds=10, seed=42)
        best_nrounds = cv_results.shape[0] - 1
        
        # 用最佳迭代次数训练模型
        bst = xgb.train(params, dtrain, num_boost_round=best_nrounds)
        dtest = xgb.DMatrix(X_test)
        y_pred_cv = bst.predict(dtest)
        
        # 计算RMSE
        rmse_cv = np.sqrt(metrics.mean_squared_error(y_test, y_pred_cv))
        
        # 打印当前参数组合和RMSE
        print(f"Params after CV: {params}, RMSE after CV: {rmse_cv}")
        
        # 记录最优参数组合
        if rmse_cv < best_rmse:
            best_rmse = rmse_cv
            best_params = params

print(f"Best params after CV: {best_params}, Best RMSE after CV: {best_rmse}")

# 用最佳参数组合重新训练模型
# 转换数据为DMatrix格式
dtrain = xgb.DMatrix(X_train, label=y_train)
cv_results = xgb.cv(best_params, dtrain, num_boost_round=1000, nfold=10, metrics={'rmse'}, early_stopping_rounds=10, seed=42)
best_nrounds = cv_results.shape[0] - 1

# 用最佳迭代次数训练模型
bst = xgb.train(best_params, dtrain, num_boost_round=best_nrounds)
dtest = xgb.DMatrix(X_test)
y_pred_cv = bst.predict(dtest)


# 计算交叉验证后评估指标
mse_cv = metrics.mean_squared_error(y_test, y_pred_cv)
rmse_cv = np.sqrt(mse_cv)
mae_cv = metrics.mean_absolute_error(y_test, y_pred_cv)
r2_cv = metrics.r2_score(y_test, y_pred_cv)

print("交叉验证后均方误差 (MSE):", mse_cv)
print("交叉验证后均方根误差 (RMSE):", rmse_cv)
print("交叉验证后平均绝对误差 (MAE):", mae_cv)
print("交叉验证后拟合优度 (R-squared):", r2_cv)

# 可视化特征重要性
plt.figure()
plot_importance(bst, ax=plt.gca())
plt.title('Feature Importance')
plt.xlabel('F-Score')  # 设置横坐标轴标题
plt.ylabel('Features')  # 设置纵坐标轴标题
plt.savefig('feature_importance_cv_TPM_1.0.png')  # 保存特征重要性图
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

plt.savefig('scatter_plot_TPM_1.0.png')  # 保存散点图
# 显示图表
plt.show()

# 保存模型
bst.save_model('my_model_TPM_1.0.json')
