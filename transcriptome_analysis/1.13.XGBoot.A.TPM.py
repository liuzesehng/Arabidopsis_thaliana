import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn.model_selection import GridSearchCV
from sklearn import metrics
from xgboost import plot_importance
import joblib
from sklearn.impute import KNNImputer
import shap
import os

# 创建目录（如果不存在）
output_dir = "TPM_1.4"
os.makedirs(output_dir, exist_ok=True)
#读取数据
df = pd.read_csv('/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA3/Alt.RCA.tsv', sep='\t')

# 检查是否还有缺失值
df.isnull().sum()

# 去除包含缺失值的行
#data = df.dropna()

# 划分特征和目标变量
x = df.drop(['TPM', 'α_TPM/%', 'β_TPM/%', 'β1_TPM/%', 'β2_TPM/%'], axis=1)
x = x.apply(pd.to_numeric, errors='coerce')
y = df['α_TPM/%']

# 检查并处理缺失值
x = x.fillna(x.mean())
y = y.fillna(y.mean())

# 划分训练集和测试集
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# XGBoost模型参数
params_xgb = {
    'learning_rate': 0.02,            # 学习率
    'booster': 'gbtree',              # 提升方法，这里使用梯度提升树（Gradient Boosting Tree）
    'objective': 'reg:squarederror',  # 损失函数，这里使用平方误差 
    'max_leaves': 127,                # 每棵树的叶子节点数量，控制模型复杂度。较大值可以提高模型复杂度但可能导致过拟合    
    'verbosity': 1,                   # 控制 XGBoost 输出信息的详细程度，0表示无输出，1表示输出进度信息    
    'seed': 42,                       # 随机种子，用于重现模型的结果    
    'nthread': 8,                    # 并行运算的线程数量，-1表示使用所有可用的CPU核心    
    'colsample_bytree': 0.6,          # 每棵树随机选择的特征比例，用于增加模型的泛化能力    
    'subsample': 0.7                  # 每次迭代时随机选择的样本比例，用于增加模型的泛化能力
}

#初始化回归模型
model_xgb = xgb.XGBRegressor(**params_xgb)

# 生成 max_leaves 的参数范围列表，从 10 到 500
max_leaves_range = list(range(120, 130))

# 定义参数网格，用于网格搜索
param_grid = {    
    'n_estimators': [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000],  # 树的数量，控制模型的复杂度    
    'max_depth': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],  # 树的最大深度，控制模型的复杂度，防止过拟合    
    'min_child_weight': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],  # 节点最小权重，值越大，算法越保守，用于控制过拟合
    'learning_rate': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],   # 学习率    
    'subsample': [0.6, 0.7, 0.8, 0.9],  # 每次迭代时随机选择的样本比例，防止过拟合    
    'colsample_bytree': [0.5, 0.6, 0.7, 0.8],  # 每棵树随机选择的特征比例，防止过拟合
    'max_leaves': max_leaves_range,     # 每棵树的叶子节点数量，控制模型复杂度。较大值可以提高模型复杂度但可能导致过拟合
}

#使用GridSearchCV进行网格搜索和k折交叉验证
grid_search = GridSearchCV(    
    estimator=model_xgb,    
    param_grid=param_grid,    
    scoring='neg_root_mean_squared_error',    
    cv=5,    
    n_jobs=8,    
    verbose=1
)

# 训练模型
grid_search.fit(x_train, y_train)

# 输出最优参数
print("Best parameters found: ", grid_search.best_params_)
print("Best RMSE score: ", -grid_search.best_score_)

# 使用最优参数训练模型
best_model = grid_search.best_estimator_

# 预测
y_pred = best_model.predict(x_test)
y_pred_list = y_pred.tolist()  
mse = metrics.mean_squared_error(y_test, y_pred_list)
rmse = np.sqrt(mse)
mae = metrics.mean_absolute_error(y_test, y_pred_list)
r2 = metrics.r2_score(y_test, y_pred_list)
print("均方误差 (MSE):", mse)
print("均方根误差 (RMSE):", rmse)
print("平均绝对误差 (MAE):", mae)
print("拟合优度 (R-squared):", r2)

# 计算shap值
explainer = shap.TreeExplainer(best_model)
# 计算shap值为numpy.array数组
shap_values_numpy = explainer.shap_values(x)
shap_values_numpy

# 计算shap值为Explanation格式
shap_values_Explanation = explainer(x)

# 创建主图（用来画蜂巢图）
fig, ax1 = plt.subplots(figsize=(10, 8), dpi=1200)
# 在主图上绘制蜂巢图，并保留热度条
shap.summary_plot(shap_values_numpy, x, feature_names=x.columns, plot_type="dot", show=False, color_bar=True)
plt.gca().set_position([0.2, 0.2, 0.65, 0.65])  #调整图表位置，留出右侧空间放热度条
# 获取共享的 y 轴
ax1 = plt.gca()
# 创建共享 y 轴的另一个图，绘制特征贡献图在顶部x轴
ax2 = ax1.twiny()
shap.summary_plot(shap_values_numpy, x, plot_type="bar", show=False)
plt.gca().set_position([0.2, 0.2, 0.65, 0.65])  # 调整图表位置，与蜂巢图对齐
# 在顶部 X 轴添加一条横线
ax2.axhline(y=13, color='gray', linestyle='-', linewidth=1)  # 注意y值应该对应顶部
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
plt.savefig(f"{output_dir}/SHAP_corrected_A_TPM_1.0.pdf", format='pdf', bbox_inches='tight')
plt.show()

# 绘制 'Tem' 特征的SHAP依赖图
shap.dependence_plot('Tem/℃', shap_values_Explanation.values, x, interaction_index='Lat', show=False)
plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)
plt.savefig(f"{output_dir}/SHAP_Dependence_A_TPM_1.0.pdf", format='pdf', bbox_inches='tight', dpi=1200)

# 绘制 'Lat' 特征的SHAP依赖图
shap.dependence_plot('Lat', shap_values_Explanation.values, x, interaction_index='Tem/℃', show=False)
plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)
plt.savefig(f"{output_dir}/SHAP_Dependence_A_Lat_TPM_1.0.pdf", format='pdf', bbox_inches='tight', dpi=1200)

# 绘制 'Lon' 特征的SHAP依赖图
shap.dependence_plot('Long', shap_values_Explanation.values, x, show=False)
plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)
plt.savefig(f"{output_dir}/SHAP_Dependence_A_Long_TPM_1.0.pdf", format='pdf', bbox_inches='tight', dpi=1200)

# 绘制 'Alt' 特征的SHAP依赖图
shap.dependence_plot('alt/m', shap_values_Explanation.values, x, show=False)
plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)
plt.savefig(f"{output_dir}/SHAP_Dependence_A_Alt_TPM_1.0.pdf", format='pdf', bbox_inches='tight', dpi=1200)

# 可视化特征重要性
plt.figure()
plot_importance(best_model, ax=plt.gca())
plt.title('Feature Importance')
plt.xlabel('F-Score')  # 设置横坐标轴标题
plt.ylabel('Features')  # 设置纵坐标轴标题
plt.savefig(f'{output_dir}/feature_A_TPM_1.0.pdf', format='pdf', bbox_inches='tight', dpi=1200)  # 保存特征重要性图
plt.show()

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

plt.title('AT rubisco activase A')
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.legend(loc="upper left")
#plt.grid(True) # 添加网格
plt.savefig(f'{output_dir}/A_TPM_1.0.pdf', format='pdf', bbox_inches='tight', dpi=1200)  # 保存特征重要性图
plt.show()


# 保存模型
best_model.save_model(f'{output_dir}/my_model_A_TPM_1.0.json')
joblib.dump(best_model , f'{output_dir}/XGBoost.A_TPM_1.0.pkl')


