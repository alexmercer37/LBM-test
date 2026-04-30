import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.interpolate import griddata
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
# 2. 读取数据
# ===============================
def load_data(file_path):
    df = pd.read_csv(file_path)
    df = df.dropna(subset=['x', 'y', 'ux', 'uy', 'phi'])
    df['x_norm'] = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    df['u_mag'] = np.sqrt(df['ux']**2 + df['uy']**2)
    df['u_mag_norm'] = (df['u_mag'] - df['u_mag'].min()) / (df['u_mag'].max() - df['u_mag'].min())
    return df

# ===============================
# 3. 绘制单帧流场（横线密集，竖线只出现在紧贴壁面的极薄层）
# ===============================
def draw_frame(df, step_num, output_path, 
               grid_density=35,            # 插值网格分辨率（影响流线平滑度）
               stream_density=1.0,         # 流线密度（值越大横线越多）
               velocity_threshold=0.02,    # 低速阈值
               wall_uy_thickness=0.02):    # 保留 uy 的壁面薄层厚度（径向比例）
    fig, ax = plt.subplots(figsize=(12, 5))
    
    x_norm = df['x_norm'].values
    y = df['y'].values
    u_mag_norm = df['u_mag_norm'].values
    ux = df['ux'].values
    uy = df['uy'].values
    phi = df['phi'].values
    
    # 速度幅值云图
    triang = tri.Triangulation(x_norm, y)
    tcf = ax.tripcolor(triang, u_mag_norm, shading='gouraud', cmap='jet', vmin=0, vmax=1)
    cbar = fig.colorbar(tcf, ax=ax, label='Velocity magnitude')
    cbar.set_label('Velocity magnitude', fontsize=12)
    
    # 创建规则网格
    xi = np.linspace(0, 1, grid_density)
    yi = np.linspace(y.min(), y.max(), grid_density)
    Xg, Yg = np.meshgrid(xi, yi)
    
    # 插值速度分量
    ux_grid = griddata((x_norm, y), ux, (Xg, Yg), method='linear')
    uy_grid = griddata((x_norm, y), uy, (Xg, Yg), method='linear')
    
    ux_grid = np.nan_to_num(ux_grid, nan=0.0)
    uy_grid = np.nan_to_num(uy_grid, nan=0.0)
    
    # ========== 关键：只在紧贴壁面的极薄层内保留径向速度 ==========
    y_range = y.max() - y.min()
    y_low  = y.min() + wall_uy_thickness * y_range
    y_high = y.max() - wall_uy_thickness * y_range
    # 将中间大部分区域的 uy 清零，只保留壁面附近薄层内的 uy
    uy_grid[(Yg >= y_low) & (Yg <= y_high)] = 0.0
    
    # 低速阈值处理（清除极低速度，避免杂线）
    u_mag_grid = np.sqrt(ux_grid**2 + uy_grid**2)
    max_u = u_mag_grid.max()
    if max_u > 0:
        low_mask = u_mag_grid < (velocity_threshold * max_u)
        ux_grid[low_mask] = 0.0
        uy_grid[low_mask] = 0.0
    
    # 绘制流线（横线会非常密集）
    ax.streamplot(Xg, Yg, ux_grid, uy_grid,
                  color='white', linewidth=0.5,
                  density=stream_density,
                  arrowstyle='->', arrowsize=0.8,
                  broken_streamlines=True)
    
    # 相区分界线（可选）
    if phi.max() >= 0.90:
        ax.tricontour(triang, phi, levels=[0.90], colors='white', linewidths=1, linestyles='--', alpha=0.7)
    if phi.min() <= 0.02:
        ax.tricontour(triang, phi, levels=[0.02], colors='white', linewidths=1, linestyles='--', alpha=0.7)
    
    # 简短的文字标注
    y_center = (y.min() + y.max()) / 2
    liquid_mask = phi > 0.90
    if np.any(liquid_mask):
        ax.text(x_norm[liquid_mask].mean(), y_center, 'Liquid', ha='center', fontsize=10, color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.6))
    gas_mask = phi < 0.02
    if np.any(gas_mask):
        ax.text(x_norm[gas_mask].mean(), y_center, 'Gas', ha='center', fontsize=10, color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.6))
    two_mask = (phi >= 0.02) & (phi <= 0.90)
    if np.any(two_mask):
        if np.any(liquid_mask) and np.any(gas_mask):
            x_tw = (x_norm[liquid_mask].mean() + x_norm[gas_mask].mean())/2
        else:
            x_tw = x_norm[two_mask].mean()
        ax.text(x_tw, y_center, 'Two-Phase', ha='center', fontsize=10, color='white',
                bbox=dict(boxstyle='round', fc='black', alpha=0.6))
    
    ax.set_title(f'Velocity Field (Step {step_num})', fontsize=14)
    ax.set_xlabel('Axial position X', fontsize=12)
    ax.set_ylabel('Radial position Y', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(y.min(), y.max())
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

# ===============================
# 4. 主程序
# ===============================
def main():
    data_dir = './'
    temp_dir = './velocity_streamlines'
    os.makedirs(temp_dir, exist_ok=True)
    
    files = get_sorted_files(data_dir, 'results_*.csv')
    if not files:
        print("未找到 results_*.csv 文件！")
        return
    
    print(f"找到 {len(files)} 个文件，生成流线动画（横线密集，竖线仅紧贴壁面）...")
    frame_files = []
    for idx, f in enumerate(files):
        step_num = re.findall(r'\d+', os.path.basename(f))[0]
        df = load_data(f)
        out_path = os.path.join(temp_dir, f'frame_{idx:04d}.png')
        # 可根据需要调整参数
        draw_frame(df, step_num, out_path,
                   grid_density=35,
                   stream_density=1.0,
                   velocity_threshold=0.02,
                   wall_uy_thickness=0.02)   # 竖线只出现在上下各2%的壁面薄层内
        frame_files.append(out_path)
        if (idx+1) % 10 == 0:
            print(f"已生成 {idx+1}/{len(files)} 帧")
    
    print("合成 GIF...")
    images = [imageio.imread(f) for f in frame_files]
    imageio.mimsave('velocity_streamlines.gif', images, fps=5, loop=0)
    print("保存为 velocity_streamlines.gif")

if __name__ == "__main__":
    main()