import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import glob
import re
import os
import imageio

# ===============================
# 1. 获取所有 CSV 文件并排序
# ===============================
def get_sorted_files(data_dir='./', pattern='results_*.csv'):
    files = glob.glob(os.path.join(data_dir, pattern))
    files.sort(key=lambda f: int(re.findall(r'\d+', os.path.basename(f))[0]))
    return files

# ===============================
# 2. 读取单个 CSV 并返回处理后的 DataFrame（保留原始温度 T）
# ===============================
def load_data(file_path):
    df = pd.read_csv(file_path, usecols=range(5))  # x, y, rho, T, phi
    df = df.dropna(subset=['x', 'y', 'T', 'phi'])
    # 归一化 x 坐标到 0~1
    df['x_norm'] = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    # 不再归一化温度，直接保留原始 T
    return df

# ===============================
# 3. 为单个文件绘制静态图片（温度场，使用全局统一的颜色范围）
# ===============================
def draw_frame(df, step_num, output_path, vmin, vmax):
    fig, ax = plt.subplots(figsize=(12, 5))
    
    x_norm = df['x_norm'].values
    y = df['y'].values
    T = df['T'].values           # 原始温度
    phi = df['phi'].values
    
    triang = tri.Triangulation(x_norm, y)
    
    # 温度云图，使用固定的 vmin, vmax
    tcf = ax.tripcolor(triang, T, shading='gouraud', cmap='jet', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(tcf, ax=ax, label='Temperature')
    cbar.set_label('Temperature', fontsize=12)
    
    # 相区分界线
    if phi.max() >= 0.90:
        ax.tricontour(triang, phi, levels=[0.90],
                      colors='white', linewidths=1.5, linestyles='--')
    if phi.min() <= 0.02:
        ax.tricontour(triang, phi, levels=[0.02],
                      colors='white', linewidths=1.5, linestyles='--')
    
    # 动态文字
    y_center = (y.min() + y.max()) / 2
    
    liquid_mask = phi > 0.90
    if np.any(liquid_mask):
        liquid_x_center = x_norm[liquid_mask].mean()
        ax.text(liquid_x_center, y_center, 'Liquid\nRegion',
                ha='center', va='center', fontsize=11, color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.7))
    
    gas_mask = phi < 0.02
    if np.any(gas_mask):
        gas_x_center = x_norm[gas_mask].mean()
        ax.text(gas_x_center, y_center, 'Gas\nRegion',
                ha='center', va='center', fontsize=11, color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.7))
    
    two_phase_mask = (phi >= 0.02) & (phi <= 0.90)
    if np.any(two_phase_mask):
        if np.any(liquid_mask) and np.any(gas_mask):
            two_phase_x_center = (x_norm[liquid_mask].mean() + x_norm[gas_mask].mean()) / 2
        else:
            two_phase_x_center = x_norm[two_phase_mask].mean()
        ax.text(two_phase_x_center, y_center, 'Two-Phase\nRegion',
                ha='center', va='center', fontsize=11, color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.7))
    
    # 流动方向箭头
    arrow_y = y.max() - (y.max() - y.min()) * 0.1
    ax.annotate('', xy=(0.9, arrow_y), xytext=(0.1, arrow_y),
                arrowprops=dict(arrowstyle='->', color='white', lw=2))
    ax.text(0.5, arrow_y + (y.max() - y.min()) * 0.02,
            'Flow Direction', ha='center', fontsize=10, color='white')
    
    # 坐标轴
    ax.set_title(f'Temperature Field (Step {step_num})', fontsize=14)
    ax.set_xlabel('Axial position X', fontsize=12)
    ax.set_ylabel('Radial position Y', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(y.min(), y.max())
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

# ===============================
# 4. 主程序：先计算全局温度范围，再逐帧绘图并合成动画
# ===============================
def main():
    data_dir = './'                # CSV 文件目录
    temp_dir = './temp_frames'     # 临时图片目录
    os.makedirs(temp_dir, exist_ok=True)
    
    files = get_sorted_files(data_dir, 'results_*.csv')
    if not files:
        print("未找到任何 results_*.csv 文件！")
        return
    
    # 第一步：遍历所有文件，获取全局温度最小值与最大值
    print("正在扫描所有文件以获取全局温度范围...")
    T_min = float('inf')
    T_max = -float('inf')
    for file_path in files:
        df = pd.read_csv(file_path, usecols=range(5))
        T_min = min(T_min, df['T'].min())
        T_max = max(T_max, df['T'].max())
    print(f"全局温度范围: {T_min:.3f} ~ {T_max:.3f}")
    
    # 可选：根据物理常识微调范围（例如入口 1.7，壁面 2.2，饱和 0.5）
    # 这里直接使用数据中的极值，但也可以手动设置为 T_min=0.5, T_max=2.2
    # 如果希望严格按给定的边界条件，取消下面两行注释：
    # T_min = 0.5
    # T_max = 2.2
    
    print(f"找到 {len(files)} 个文件，开始生成温度场图片帧...")
    
    frame_files = []
    for idx, file_path in enumerate(files):
        step_num = re.findall(r'\d+', os.path.basename(file_path))[0]
        df = load_data(file_path)
        out_path = os.path.join(temp_dir, f'frame_{idx:04d}.png')
        draw_frame(df, step_num, out_path, vmin=T_min, vmax=T_max)
        frame_files.append(out_path)
        if (idx+1) % 10 == 0:
            print(f"已生成 {idx+1}/{len(files)} 帧")
    
    print("所有帧生成完毕，开始合成 GIF...")
    
    images = []
    for fname in frame_files:
        images.append(imageio.imread(fname))
    
    gif_path = 'temperature_evolution.gif'
    imageio.mimsave(gif_path, images, fps=5, loop=0)
    print(f"温度场动画已保存为 {gif_path}")

if __name__ == "__main__":
    main()