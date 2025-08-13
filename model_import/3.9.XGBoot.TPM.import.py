import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn import metrics
import joblib
import shap
import os
from sklearn.feature_selection import SelectFromModel
from xgboost import XGBRegressor

# 创建输出目录（如果不存在）
output_dir = "TPM_import_results"
os.makedirs(output_dir, exist_ok=True)

# 加载已训练好的模型
print("正在加载预训练模型...")
try:
    # 方法1：使用XGBoost原生方法加载模型
    model_path = "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_3.4/my_model_TPM_3.0.json"
    loaded_model = xgb.XGBRegressor()
    loaded_model.load_model(model_path)
    print(f"成功从 {model_path} 加载模型")
except Exception as e:
    print(f"使用JSON格式加载模型失败: {e}")
    try:
        # 方法2：使用joblib加载模型（如果JSON格式失败）
        pkl_model_path = "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_3.4/TPM_XGBoost_3.0.pkl"
        loaded_model = joblib.load(pkl_model_path)
        print(f"成功从 {pkl_model_path} 加载模型")
    except Exception as e2:
        print(f"使用PKL格式加载模型也失败: {e2}")
        exit(1)

# 读取数据（使用与训练时相同的数据）
print("正在读取数据...")
df = pd.read_csv('/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA3/Alt.RCA.me_snp.tsv', sep='\t')

# 检查数据
print(f"数据形状: {df.shape}")
print(f"数据列数: {len(df.columns)}")

# 划分特征和目标变量（与训练时保持一致）
x = df.drop(['TPM', 'α_TPM/%', 'β_TPM/%', 'β1_TPM/%', 'β2_TPM/%'], axis=1)
x = x.apply(pd.to_numeric, errors='coerce')
y = df['TPM']

# 处理缺失值（与训练时保持一致）
x = x.fillna(x.mean())
y = y.fillna(y.mean())

# 划分训练集和测试集
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)

# 读取训练时选择的特征（如果存在）
selected_features_file = "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_3.4/selected_features_TPM_3.0.txt"

if os.path.exists(selected_features_file):
    print("正在加载训练时选择的特征...")
    with open(selected_features_file, 'r') as f:
        selected_features_str = f.read().strip()
        selected_features = selected_features_str.split('\t')
    
    # 使用训练时选择的特征
    X_train_selected = x_train[selected_features]
    X_test_selected = x_test[selected_features]
    print(f"使用训练时选择的 {len(selected_features)} 个特征")
else:
    print("未找到特征选择文件，使用所有特征进行预测...")
    X_train_selected = x_train
    X_test_selected = x_test
    selected_features = x.columns.tolist()

# 使用加载的模型进行预测
print("正在进行预测...")
try:
    y_pred = loaded_model.predict(X_test_selected)
    y_pred_list = y_pred.tolist()  
    print("预测完成!")
    
    # 计算评估指标
    mse = metrics.mean_squared_error(y_test, y_pred_list)
    rmse = np.sqrt(mse)
    mae = metrics.mean_absolute_error(y_test, y_pred_list)
    r2 = metrics.r2_score(y_test, y_pred_list)
    
    print("均方误差 (MSE):", mse)
    print("均方根误差 (RMSE):", rmse)
    print("平均绝对误差 (MAE):", mae)
    print("拟合优度 (R-squared):", r2)
    
    
    # 计算SHAP值进行特征重要性分析
    try:
        print("\n正在计算SHAP值...")
        explainer = shap.TreeExplainer(loaded_model)

        # 从完整数据集 x 中选择相同的特征
        X_selected = x[selected_features]

        # 计算shap值为numpy.array数组
        shap_values_numpy = explainer.shap_values(X_selected)

        # 计算shap值为Explanation格式
        shap_values_Explanation = explainer(X_selected)

        print("SHAP分析完成!")

        # 遍历不同类型的甲基化特征
        feature_types = ["CG_promoter", "CG", "CG_terminator", "CHG_promoter", "CHG", "CHG_terminator", "CHH_promoter", "CHH", "CHH_terminator"]
        
        for feature_type in feature_types:
            if feature_type in selected_features:
                # 绘制 SHAP 依赖图
                shap.dependence_plot(feature_type, shap_values_Explanation.values, X_selected, 
                                   feature_names=selected_features, show=False)    
                # 添加水平线
                plt.axhline(y=0, color='black', linestyle='-.', linewidth=1)   
                # 保存图像
                plt.savefig(f"{output_dir}/SHAP_Dependence_{feature_type}_TPM_3.0.pdf", format='pdf', bbox_inches='tight', dpi=1200)
                # 清除当前图像，以便绘制下一个特征
                plt.clf()
                
    except Exception as shap_error:
        print(f"SHAP分析失败: {shap_error}")
        print("跳过SHAP分析，继续其他分析...")    

except Exception as pred_error:
    print(f"预测过程中出现错误: {pred_error}")
    print("请检查模型文件和数据文件是否匹配")

print("\n模型导入和预测分析完成!")
