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
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

BASE_FONT_SIZE = 16
LABEL_FONT_SIZE = 18
TITLE_FONT_SIZE = 20
TICK_FONT_SIZE = 15
LEGEND_FONT_SIZE = 15

plt.rcParams.update({
    'font.size': BASE_FONT_SIZE,
    'axes.titlesize': TITLE_FONT_SIZE,
    'axes.labelsize': LABEL_FONT_SIZE,
    'xtick.labelsize': TICK_FONT_SIZE,
    'ytick.labelsize': TICK_FONT_SIZE,
    'legend.fontsize': LEGEND_FONT_SIZE
})

# 定义变量名
name_TPM = "A_TPM"

# 读取原始数据
df = pd.read_csv('../RCA/Alt.snp_meth.filtered.tsv', sep='\t')

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
x = df.drop(['Temperature', 'total', 'α', 'β', 'β1', 'β2', "α/%", "β/%", "β1/%", "β2/%"], axis=1)
x = x.apply(pd.to_numeric, errors='coerce')
y = df['α']

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

# 计算测试集shap值
explainer = shap.TreeExplainer(best_model)
shap_values_numpy = explainer.shap_values(x_test)
shap_values_Explanation = explainer(x_test)
    
# 获取前10个最重要的特征
mean_abs_shap = np.abs(shap_values_numpy).mean(axis=0)
top_10_indices = np.argsort(mean_abs_shap)[-10:][::-1]
top_10_features = x_test.columns[top_10_indices].tolist()

threshold_value = 0.1

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
threshold_indices = sorted_indices[mean_abs_shap[sorted_indices] > threshold_value]
threshold_features = x_test.columns[threshold_indices].tolist()
threshold_feature_values = mean_abs_shap[threshold_indices]
sorted_contributions = mean_abs_shap[sorted_indices]
cumulative_contribution = np.cumsum(sorted_contributions) / total_contribution * 100

print(f"\nmean(|SHAP|) > {threshold_value} 的特征数量: {len(threshold_features)}")
print(f"特征列表: {threshold_features}")

# 找到累计贡献度达到80%的特征索引
threshold_80_idx = np.where(cumulative_contribution >= 80)[0][0]
top_80_features = x_test.columns[sorted_indices[:threshold_80_idx + 1]].tolist()
top_80_contribution = cumulative_contribution[threshold_80_idx]

print(f"\n累计贡献度达到80%的特征数量: {len(top_80_features)}")
print(f"累计贡献度: {top_80_contribution:.2f}%")
print(f"特征列表: {top_80_features}")

# 将阈值特征信息追加到文件
with open(f'{output_dir}/{name_TPM}_optimization_results.txt', 'a', encoding='utf-8') as f:
    f.write(f"\nSHAP阈值特征分析:\n")
    f.write(f"  说明: 基于测试集 SHAP，阈值 = {threshold_value}\n")
    f.write(f"  mean(|SHAP|) > {threshold_value} 的特征数量: {len(threshold_features)}\n")
    f.write(f"  特征列表: {', '.join(threshold_features) if threshold_features else '无'}\n")
    f.write("  特征及其平均|SHAP|值:\n")
    if threshold_features:
        for feature, value in zip(threshold_features, threshold_feature_values):
            f.write(f"    {feature}: {value:.6f}\n")
    else:
        f.write("    无\n")

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


def beeswarm_offsets(values, nbins=60, spread=0.8):
    values = np.asarray(values)
    if len(values) <= 1:
        return np.zeros(len(values))

    bins = np.linspace(values.min(), values.max(), nbins + 1)
    bin_ids = np.digitize(values, bins) - 1
    offsets = np.zeros(len(values), dtype=float)

    for bin_id in np.unique(bin_ids):
        idx = np.where(bin_ids == bin_id)[0]
        if len(idx) <= 1:
            continue

        sorted_idx = idx[np.argsort(values[idx], kind='mergesort')]
        pattern = np.arange(len(sorted_idx), dtype=float)
        pattern = np.where(pattern % 2 == 0, pattern / 2, -(pattern + 1) / 2)
        max_abs = np.max(np.abs(pattern))
        if max_abs > 0:
            pattern = pattern / max_abs * spread
        offsets[sorted_idx] = pattern

    return offsets


max_display = min(10, x_test.shape[1])
plot_indices = sorted_indices[:max_display]
plot_features = x_test.columns[plot_indices].tolist()
plot_mean_abs_shap = mean_abs_shap[plot_indices]

threshold_plot_count = int(np.sum(plot_mean_abs_shap > threshold_value))
threshold_line_y = threshold_plot_count - 0.5 if 0 < threshold_plot_count < max_display else None

cmap = plt.get_cmap('cool')
fig = plt.figure(figsize=(12, 6), dpi=1200)
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.9, 0.05], wspace=0.15)
ax_bar = fig.add_subplot(gs[0, 0])
ax_bee = fig.add_subplot(gs[0, 1], sharey=ax_bar)
cax = fig.add_subplot(gs[0, 2])
y_positions = np.arange(max_display)

ax_bar.barh(y_positions, plot_mean_abs_shap, color='#2f83d0')
ax_bar.set_yticks(y_positions)
ax_bar.set_yticklabels(plot_features, fontsize=TICK_FONT_SIZE)
ax_bar.invert_yaxis()
ax_bar.set_xlabel('Mean(|SHAP value|)', fontsize=LABEL_FONT_SIZE)
ax_bar.grid(axis='y', linestyle=':', linewidth=0.8, alpha=0.4)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.tick_params(axis='x', labelsize=TICK_FONT_SIZE)
ax_bar.tick_params(axis='y', length=0, labelsize=TICK_FONT_SIZE)

for row_idx, feature_idx in enumerate(plot_indices):
    feature_name = x_test.columns[feature_idx]
    feature_values = x_test[feature_name].to_numpy()
    shap_row_values = shap_values_numpy[:, feature_idx]
    offsets = beeswarm_offsets(shap_row_values)

    vmin = np.nanpercentile(feature_values, 5)
    vmax = np.nanpercentile(feature_values, 95)
    if np.isclose(vmin, vmax):
        vmin = np.nanmin(feature_values)
        vmax = np.nanmax(feature_values)
    if np.isclose(vmin, vmax):
        vmin -= 1
        vmax += 1
    norm = Normalize(vmin=vmin, vmax=vmax)

    ax_bee.scatter(
        shap_row_values,
        np.full_like(shap_row_values, row_idx, dtype=float) + offsets * 0.32,
        c=feature_values,
        cmap=cmap,
        norm=norm,
        s=16,
        alpha=0.95,
        linewidths=0
    )

ax_bee.axvline(0, color='gray', linewidth=1)
ax_bee.set_yticks(y_positions)
ax_bee.tick_params(axis='x', labelsize=TICK_FONT_SIZE)
ax_bee.tick_params(axis='y', left=False, labelleft=False, labelsize=TICK_FONT_SIZE)
ax_bee.set_xlabel('SHAP value (impact on model output)', fontsize=LABEL_FONT_SIZE)
ax_bee.grid(axis='y', linestyle=':', linewidth=0.8, alpha=0.4)
ax_bee.spines['top'].set_visible(False)
ax_bee.spines['right'].set_visible(False)
ax_bee.spines['left'].set_visible(False)

if threshold_line_y is not None:
    ax_bar.axhline(y=threshold_line_y, color='black', linewidth=1.1)
    ax_bee.axhline(y=threshold_line_y, color='black', linewidth=1.1)

colorbar = fig.colorbar(ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap=cmap), cax=cax)
colorbar.set_ticks([0, 1])
colorbar.set_ticklabels(['Low', 'High'])
colorbar.ax.tick_params(labelsize=TICK_FONT_SIZE)
colorbar.set_label('Feature value', rotation=270, labelpad=20, fontsize=LABEL_FONT_SIZE)
colorbar.outline.set_visible(False)

fig.savefig(f"{output_dir}/SHAP_combined_{name_TPM}_4.0.pdf",
            format='pdf', bbox_inches='tight')
plt.close(fig)

# 绘制前10个最重要特征的SHAP依赖图
for feature in top_10_features:
    # 绘制依赖图，interaction_index默认为'auto'，会自动寻找交互最强的特征
    shap.dependence_plot(
    feature,
    shap_values_Explanation.values,
    x_test,
    interaction_index='auto',
    x_jitter=0,
    dot_size=12,
    alpha=0.9,
    show=False
)
    ax = plt.gca()
    feature_min = x_test[feature].min()
    feature_max = x_test[feature].max()
    if feature_min >= 0:
        right_margin = (feature_max - feature_min) * 0.03 if feature_max > feature_min else 0.1
        ax.set_xlim(left=0, right=feature_max + right_margin)
    ax.set_xlabel(ax.get_xlabel(), fontsize=LABEL_FONT_SIZE)
    ax.set_ylabel(ax.get_ylabel(), fontsize=LABEL_FONT_SIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONT_SIZE)
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

plt.title(r'Expression level of Rca $\alpha$', fontsize=TITLE_FONT_SIZE)
plt.xlabel('Actual Values', fontsize=LABEL_FONT_SIZE)
plt.ylabel('Predicted Values', fontsize=LABEL_FONT_SIZE)
plt.xticks(fontsize=TICK_FONT_SIZE)
plt.yticks(fontsize=TICK_FONT_SIZE)
plt.legend(loc="upper left", fontsize=LEGEND_FONT_SIZE)
#plt.grid(True) # 添加网格
plt.savefig(f'{output_dir}/{name_TPM}_4.0.pdf', format='pdf', bbox_inches='tight', dpi=1200)  # 保存特征重要性图
plt.show()
plt.close()


# 保存模型
best_model.save_model(f'{output_dir}/my_model_{name_TPM}_4.0.json')
