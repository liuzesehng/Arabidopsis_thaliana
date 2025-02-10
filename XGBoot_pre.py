import pandas as pd
import numpy as np
import xgboost as xgb
import shap

# 加载保存的模型
model = xgb.Booster()
model.load_model('my_model_TPM_1.0.json')

# 数据读取
df = pd.read_csv('../RCA/RCA.tsv', sep='\t')
df.head()

# 缺失值处理
df.isnull().sum()
# 划分特征和目标变量
explainer = shap.TreeExplainer(model)
# 计算shap值为numpy.array数组
shap_values_numpy = explainer.shap_values(x)
shap_values_numpy


# 假设这是新的数据样本（一个或多个样本），需要与训练数据的特征相同
new_data = np.array([[8.3252, 41.0, 6.984127, 1.023810, 322.0, 2.555556, 37.88, -122.23],
                     [8.3014, 21.0, 6.238137, 0.971880, 2401.0, 2.109842, 37.86, -122.22]])

# 将新数据转换为DMatrix对象，这是XGBoost使用的数据格式
dnew = xgb.DMatrix(new_data)

# 使用模型进行预测
new_predictions = model.predict(dnew)

print("New predictions:", new_predictions)