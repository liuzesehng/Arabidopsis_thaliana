import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn.model_selection import cross_val_score, KFold
from sklearn import metrics
from xgboost import plot_importance
from sklearn.impute import KNNImputer
import shap
import os
import optuna

# 定义变量名
name_TPM = "B1_per_TPM"

# 读取原始数据
df = pd.read_csv('../RCA/Alt.snp_meth.tsv', sep='\t')

# log转换
cols_to_log = ['total', 'α', 'β', 'β1', 'β2']
# 确保只转换存在的列
cols_to_transform = [col for col in cols_to_log if col in df.columns]

# 使用 log1p (log(x+1)) 避免 log(0) 导致的错误
df[cols_to_transform] = np.log1p(df[cols_to_transform])

# 创建目录（如果不存在），区分不同数据集
output_dir = f"TPM_4.5"
os.makedirs(output_dir, exist_ok=True)
   
# 检查是否还有缺失值
df.isnull().sum()

# 去除包含缺失值的行
#data = df.dropna()

# 划分特征和目标变量
x = df.drop(['tem', 'total', 'α', 'β', 'β1', 'β2', "α/%", "β/%", "β1/%", "β2/%"], axis=1)
x = x.apply(pd.to_numeric, errors='coerce')
y = df['β1/%']

# 划分训练集和测试集
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# 定义 Optuna 目标函数
def objective(trial):
    # 定义搜索空间
    param = {
        'booster': 'gbtree',
        'objective': 'reg:squarederror',
        'tree_method': 'hist',            # 加速训练
        # 'device': 'cuda',                 # 使用 GPU 加速（已注释）
        'verbosity': 0,
        'seed': 42,
        'nthread': 16,
            
        # 核心调参区间
        'n_estimators': trial.suggest_int('n_estimators', 100, 2000),
        'learning_rate': trial.suggest_float('learning_rate', 0.001, 0.2, log=True),
        'max_depth': trial.suggest_int('max_depth', 1, 15),
        'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
        'subsample': trial.suggest_float('subsample', 0.5, 1.0),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
        'reg_alpha': trial.suggest_float('reg_alpha', 1e-5, 10, log=True),
        'reg_lambda': trial.suggest_float('reg_lambda', 1e-5, 10, log=True),
    }
        
    model = xgb.XGBRegressor(**param)
        
    # 5折交叉验证，计算 R² 作为训练指标
    cv_scores = cross_val_score(model, x_train, y_train, 
                                cv=5, 
                                scoring='r2', 
                                n_jobs=1)
        
    # 返回平均 R²
    return cv_scores.mean()

# 创建优化任务
print(f"Starting Optuna optimization...")
study = optuna.create_study(direction='maximize')  # 最大化 R²
# --- 开启 Optuna 并行 ---
study.optimize(objective, n_trials=500, n_jobs=4, show_progress_bar=True)

# 输出最优参数
print(f"Best parameters found: ", study.best_params)
print(f"Best R² score: ", study.best_value)
    
# 将最优参数写入文件
with open(f'{output_dir}/{name_TPM}_optimization_results.txt', 'w', encoding='utf-8') as f:
    f.write(f"=== {name_TPM} 优化结果  ===\n\n")
    f.write(f"最优参数:\n")
    for param, value in study.best_params.items():
        f.write(f"  {param}: {value}\n")
    f.write(f"\n最优交叉验证 R² 分数: {study.best_value:.6f}\n")
    f.write(f"\n{'='*50}\n\n")

# 使用最优参数训练最终模型
best_params = {
    'booster': 'gbtree',
    'objective': 'reg:squarederror',
    'tree_method': 'hist',
    # 'device': 'cuda',  # GPU 加速（已注释）
    'verbosity': 0,
    'seed': 42,
    'nthread': 64,
    **study.best_params
}
best_model = xgb.XGBRegressor(**best_params)
best_model.fit(x_train, y_train)

# 保存 Optuna 优化历史
optuna_df = study.trials_dataframe()
optuna_df.to_csv(f'{output_dir}/{name_TPM}_optuna_optimization_history.csv', index=False)

# 预测
y_pred = best_model.predict(x_test)
y_pred_list = y_pred.tolist()  
mse = metrics.mean_squared_error(y_test, y_pred_list)
rmse = np.sqrt(mse)
mae = metrics.mean_absolute_error(y_test, y_pred_list)
r2 = metrics.r2_score(y_test, y_pred_list)
print(f"{name_TPM}预测结果:")
print("均方误差 (MSE):", mse)
print("均方根误差 (RMSE):", rmse)
print("平均绝对误差 (MAE):", mae)
print("拟合优度 (R-squared):", r2)
    
# 将预测结果追加到文件
with open(f'{output_dir}/{name_TPM}_optimization_results.txt', 'a', encoding='utf-8') as f:
    f.write(f"测试集预测结果:\n")
    f.write(f"  均方误差 (MSE): {mse:.6f}\n")
    f.write(f"  均方根误差 (RMSE): {rmse:.6f}\n")
    f.write(f"  平均绝对误差 (MAE): {mae:.6f}\n")
    f.write(f"  拟合优度 (R-squared): {r2:.6f}\n")

# 计算shap值
explainer = shap.TreeExplainer(best_model)
# 计算shap值为numpy.array数组
shap_values_numpy = explainer.shap_values(x)
    
# 计算shap值为Explanation格式
shap_values_Explanation = explainer(x)

# 创建主图（用来画蜂巢图）
fig, ax1 = plt.subplots(figsize=(10, 8), dpi=1200)
# 在主图上绘制蜂巢图，并保留热度条
shap.summary_plot(shap_values_numpy, x, feature_names=x.columns, plot_type="dot", max_display=10, show=False, color_bar=True)
plt.gca().set_position([0.2, 0.2, 0.65, 0.65])  #调整图表位置，留出右侧空间放热度条
# 获取共享的 y 轴
ax1 = plt.gca()
# 创建共享 y 轴的另一个图，绘制特征贡献图在顶部x轴
ax2 = ax1.twiny()
shap.summary_plot(shap_values_numpy, x, plot_type="bar", max_display=10, show=False)
plt.gca().set_position([0.2, 0.2, 0.65, 0.65])  # 调整图表位置，与蜂巢图对齐
# 调整透明度
bars = ax2.patches  # 获取所有的柱状图对象
for bar in bars:   
    bar.set_alpha(0.2)  # 设置透明度
# 设置两个x轴的标签
ax1.set_xlabel('Shapley Value Contribution (Bee Swarm)', fontsize=12)
ax2.set_xlabel('Mean Shapley Value (Feature Importance)', fontsize=12)
# 移动顶部的 X 轴，避免与底部 X 轴重叠
ax2.xaxis.set_label_position('top')  # 将标签移动到顶部
ax2.xaxis.tick_top()  # 将刻度也移动到顶部
# 设置y轴标签
ax1.set_ylabel('Features', fontsize=12)
plt.tight_layout()
plt.savefig(f"{output_dir}/SHAP_corrected_{name_TPM}_4.0.pdf", format='pdf', bbox_inches='tight')
plt.show()
plt.close()

# 获取前10个最重要的特征
mean_abs_shap = np.abs(shap_values_numpy).mean(axis=0)
top_10_indices = np.argsort(mean_abs_shap)[-10:][::-1]
top_10_features = x.columns[top_10_indices].tolist()

# 第一个最重要特征的贡献度
first_feature_contribution = mean_abs_shap[top_10_indices[0]]
total_contribution = mean_abs_shap.sum()
first_feature_percentage = (first_feature_contribution / total_contribution) * 100

print(f"前10个最重要的特征: {top_10_features}")
print(f"第一个最重要特征: {top_10_features[0]}")
print(f"第一个特征的平均|SHAP|值: {first_feature_contribution:.6f}")
print(f"所有特征的平均|SHAP|值总和: {total_contribution:.6f}")
print(f"第一个特征贡献度百分比: {first_feature_percentage:.2f}%")

# 计算累计贡献度达到80%的特征
sorted_indices = np.argsort(mean_abs_shap)[::-1]  # 从高到低排序
sorted_contributions = mean_abs_shap[sorted_indices]
cumulative_contribution = np.cumsum(sorted_contributions) / total_contribution * 100

# 找到累计贡献度达到80%的特征索引
threshold_80_idx = np.where(cumulative_contribution >= 80)[0][0]
top_80_features = x.columns[sorted_indices[:threshold_80_idx + 1]].tolist()
top_80_contribution = cumulative_contribution[threshold_80_idx]

print(f"\n累计贡献度达到80%的特征数量: {len(top_80_features)}")
print(f"累计贡献度: {top_80_contribution:.2f}%")
print(f"特征列表: {top_80_features}")

# 将SHAP贡献度信息追加到文件
with open(f'{output_dir}/{name_TPM}_optimization_results.txt', 'a', encoding='utf-8') as f:
    f.write(f"\nSHAP特征重要性分析:\n")
    f.write(f"  前10个最重要特征: {', '.join(top_10_features)}\n")
    f.write(f"  第一个最重要特征: {top_10_features[0]}\n")
    f.write(f"  第一个特征的平均|SHAP|值: {first_feature_contribution:.6f}\n")
    f.write(f"  所有特征的平均|SHAP|值总和: {total_contribution:.6f}\n")
    f.write(f"  第一个特征贡献度百分比: {first_feature_percentage:.2f}%\n")
    f.write(f"\n  累计贡献度达到80%的特征:\n")
    f.write(f"    特征数量: {len(top_80_features)}\n")
    f.write(f"    累计贡献度: {top_80_contribution:.2f}%\n")
    f.write(f"    特征列表: {', '.join(top_80_features)}\n")

# 绘制前10个最重要特征的SHAP依赖图
for feature in top_10_features:
    # 绘制依赖图，interaction_index默认为'auto'，会自动寻找交互最强的特征
    shap.dependence_plot(feature, shap_values_Explanation.values, x, show=False)
    plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)
        
    # 处理文件名中的特殊字符，例如将 '/' 替换为 '_'
    safe_feature_name = feature.replace('/', '_').replace('\\', '_')
    plt.savefig(f"{output_dir}/SHAP_Dependence_{safe_feature_name}_{name_TPM}_4.0.pdf", format='pdf', bbox_inches='tight', dpi=1200)
    plt.close()


# 可视化，绘制散点图
plt.figure(figsize=(8, 6), dpi=1200)
plt.scatter(y_test, y_pred, color='#f4ba8a', label="Predicted", alpha=0.2)  # 预测值散点图
# 绘制一条 y = x 的对角线，用于参考
max_value = max(max(y_test), max(y_pred))
plt.plot([0, max_value], [0, max_value], color='black', alpha=0.6, linestyle='--', label='x=y') # 1:1灰色虚线
#plt.plot(y_test, y_test, color='black', alpha=0.6, linestyle='--', label='x=y')  # 1:1灰色虚线
# 拟合线
z = np.polyfit(y_test, y_pred, 1)
p = np.poly1d(z)
plt.plot(y_test, p(y_test), color='#b4d4e1', alpha=0.6,          
            label=f"Line of Best Fit\n$R^2$ = {r2:.2f},MAE = {mae:.2f}")

plt.title(f'AT rubisco activase B1 ratio')
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.legend(loc="upper left")
#plt.grid(True) # 添加网格
plt.savefig(f'{output_dir}/{name_TPM}_4.0.pdf', format='pdf', bbox_inches='tight', dpi=1200)  # 保存特征重要性图
plt.show()
plt.close()


# 保存模型
best_model.save_model(f'{output_dir}/my_model_{name_TPM}_4.0.json')
