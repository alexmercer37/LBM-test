import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import os

# 获取所有结果文件
files = sorted(glob.glob('results_*.csv'), key=lambda f: int(re.findall(r'\d+', f)[0]))

steps = []
liquid_masses = []
vapor_masses = []
outlet_vapor_masses = []

for f in files:
    df = pd.read_csv(f)
    # 假设全局量在每一行都相同，取第一行
    # 如果 CSV 中未直接保存这些列，需要根据 phi 和 rho 重新计算（见后文备选方案）
    if 'liquid_mass' in df.columns:
        liquid_masses.append(df['liquid_mass'].iloc[0])
        vapor_masses.append(df['vapor_mass'].iloc[0])
        outlet_vapor_masses.append(df['outlet_vapor_mass'].iloc[0])
    else:
        # 备选：根据 phi 和 rho 自己计算（需要知道 rho_l, rho_g）
        # 此处仅示意，实际需从模拟参数获得
        pass
    step = int(re.findall(r'\d+', os.path.basename(f))[0])
    steps.append(step)

# 转换为 numpy 数组
steps = np.array(steps)
liquid_masses = np.array(liquid_masses)
vapor_masses = np.array(vapor_masses)
outlet_vapor_masses = np.array(outlet_vapor_masses)

# 绘制液相与气相质量演化
plt.figure(figsize=(10, 4))
plt.plot(steps, liquid_masses, label='Liquid mass', marker='o', markersize=3)
plt.plot(steps, vapor_masses, label='Vapor mass', marker='s', markersize=3)
plt.plot(steps, liquid_masses + vapor_masses, label='Total mass', linestyle='--', color='k')
plt.xlabel('Iteration step')
plt.ylabel('Mass')
plt.title('Evolution of liquid and vapor mass')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('mass_evolution.png', dpi=150)
plt.show()

# 绘制出口蒸汽质量流量
plt.figure(figsize=(10, 4))
plt.plot(steps, outlet_vapor_masses, label='Outlet vapor mass flow', color='r')
plt.xlabel('Iteration step')
plt.ylabel('Mass flow rate')
plt.title('Outlet vapor mass flow rate')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('outlet_vapor_flow.png', dpi=150)
plt.show()