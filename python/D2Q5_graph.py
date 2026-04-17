import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


filename = 'D2Q5.csv'
df = pd.read_csv(filename)

# 2. 将长格式数据转换为二维矩阵 (T 场)
# index 对应行 (i), columns 对应列 (j)
temp_matrix = df.pivot(index='x_idx', columns='y_idx', values='T')

# 3. 创建绘图
plt.figure(figsize=(12, 8))

# 使用 seaborn 绘制热力图
# cmap='magma' 或 'hot' 非常适合表现温度场
ax = sns.heatmap(temp_matrix, 
                 cmap='magma', 
                 robust=True, 
                 cbar_kws={'label': 'Temperature (T)'})

# 4. 设置标签和标题
plt.title('LBM D2Q5 Heat Conduction Simulation', fontsize=15)
plt.xlabel('Y Grid Index (n)', fontsize=12)
plt.ylabel('X Grid Index (m)', fontsize=12)


plt.show()
