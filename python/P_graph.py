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
# 2. 读取单个 CSV，返回处理后的 DataFrame
# ===============================
def load_data(file_path):
    df = pd.read_csv(file_path, usecols=range(5))  # x, y, rho, T, phi
    df = df.dropna(subset=['x', 'y', 'phi'])
    df['x_norm'] = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    return df

# ===============================
# 3. 绘制单帧相场云图
# ===============================
def draw_frame(df, step_num, output_path, vmin=0, vmax=1):
    fig, ax = plt.subplots(figsize=(12, 5))
    
    x_norm = df['x_norm'].values
    y = df['y'].values
    phi = df['phi'].values
    
    triang = tri.Triangulation(x_norm, y)
    
    # 相场云图（使用双色映射：蓝→红）
    tcf = ax.tripcolor(triang, phi, shading='gouraud', cmap='coolwarm', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(tcf, ax=ax, label='Phase field φ')
    cbar.set_label('Phase field φ', fontsize=12)
    
    # 可选：绘制相区分界线（突出液体/气体边界）
    if phi.max() >= 0.90:
        ax.tricontour(triang, phi, levels=[0.90], colors='white', linewidths=1.5, linestyles='--')
    if phi.min() <= 0.02:
        ax.tricontour(triang, phi, levels=[0.02], colors='white', linewidths=1.5, linestyles='--')
    
    # 动态文字标注（基于当前帧的 phi 分布）
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
    
    # 坐标轴与标题
    ax.set_title(f'Phase Field (Step {step_num})', fontsize=14)
    ax.set_xlabel('Axial position X', fontsize=12)
    ax.set_ylabel('Radial position Y', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(y.min(), y.max())
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

# ===============================
# 4. 主程序：逐帧保存并合成 GIF
# ===============================
def main():
    data_dir = './'               # CSV 文件所在目录
    temp_dir = './phi_frames'     # 临时图片目录
    os.makedirs(temp_dir, exist_ok=True)
    
    files = get_sorted_files(data_dir, 'results_*.csv')
    if not files:
        print("未找到任何 results_*.csv 文件！")
        return
    
    print(f"找到 {len(files)} 个文件，开始生成相场动画帧...")
    
    frame_files = []
    for idx, file_path in enumerate(files):
        step_num = re.findall(r'\d+', os.path.basename(file_path))[0]
        df = load_data(file_path)
        out_path = os.path.join(temp_dir, f'frame_{idx:04d}.png')
        # 相场值预设为 0~1，可直接使用固定 vmin=0, vmax=1
        draw_frame(df, step_num, out_path, vmin=0, vmax=1)
        frame_files.append(out_path)
        if (idx+1) % 10 == 0:
            print(f"已生成 {idx+1}/{len(files)} 帧")
    
    print("所有帧生成完毕，开始合成 GIF...")
    images = [imageio.imread(f) for f in frame_files]
    gif_path = 'phase_field_evolution.gif'
    imageio.mimsave(gif_path, images, fps=5, loop=0)
    print(f"相场动画已保存为 {gif_path}")

if __name__ == "__main__":
    main()