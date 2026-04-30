import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import glob
import re
import os
import imageio.v2 as imageio   # 使用 v2 消除警告

# ===============================
# 1. 获取所有 CSV 文件并排序
# ===============================
def get_sorted_files(data_dir='./', pattern='results_*.csv'):
    files = glob.glob(os.path.join(data_dir, pattern))
    files.sort(key=lambda f: int(re.findall(r'\d+', os.path.basename(f))[0]))
    return files

# ===============================
# 2. 读取单个 CSV，直接使用已有的 mdot 列
# ===============================
def load_data(file_path):
    df = pd.read_csv(file_path)
    required_cols = ['x', 'y', 'mdot', 'phi']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"文件 {file_path} 中缺少列 {col}")
    df = df.dropna(subset=required_cols)
    df = df.drop_duplicates(subset=['x', 'y'])
    if len(df) < 3:
        print(f"警告：{file_path} 有效点数少于3，跳过")
        return None
    df['x_norm'] = (df['x'] - df['x'].min()) / (df['x'].max() - df['x'].min())
    return df

# ===============================
# 3. 绘制单帧 mdot 云图（固定输出尺寸）
# ===============================
def draw_frame(df, step_num, output_path, vmin, vmax):
    fig, ax = plt.subplots(figsize=(12, 5))
    # 手动调整布局，确保所有帧输出尺寸一致
    fig.subplots_adjust(left=0.08, right=0.92, top=0.92, bottom=0.08)
    
    x_norm = df['x_norm'].values
    y = df['y'].values
    mdot = df['mdot'].values
    phi = df['phi'].values
    
    triang = tri.Triangulation(x_norm, y)
    
    tcf = ax.tripcolor(triang, mdot, shading='gouraud', cmap='hot', vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(tcf, ax=ax, label='Evaporation rate (mdot)')
    cbar.set_label('mdot', fontsize=12)
    
    # 相区分界线
    if phi.max() >= 0.90:
        ax.tricontour(triang, phi, levels=[0.90], colors='white', linewidths=1, linestyles='--', alpha=0.7)
    if phi.min() <= 0.02:
        ax.tricontour(triang, phi, levels=[0.02], colors='white', linewidths=1, linestyles='--', alpha=0.7)
    
    # 文字标注
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
    
    # 流动方向箭头
    # arrow_y = y.max() - (y.max() - y.min()) * 0.1
    # ax.annotate('', xy=(0.9, arrow_y), xytext=(0.1, arrow_y),
    #             arrowprops=dict(arrowstyle='->', color='white', lw=1.5))
    # ax.text(0.5, arrow_y + (y.max() - y.min()) * 0.02,
    #         'Flow', ha='center', fontsize=10, color='white')
    
    ax.set_title(f'Evaporation rate (Step {step_num})', fontsize=14)
    ax.set_xlabel('Axial position X', fontsize=12)
    ax.set_ylabel('Radial position Y', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(y.min(), y.max())
    
    # 使用 bbox_inches='tight' 可能导致尺寸微变，这里不使用 tight，直接保存
    plt.savefig(output_path, dpi=150, bbox_inches=None)
    plt.close(fig)

# ===============================
# 4. 主程序
# ===============================
def main():
    data_dir = './'
    temp_dir = './mdot_frames'
    os.makedirs(temp_dir, exist_ok=True)
    
    files = get_sorted_files(data_dir, 'results_*.csv')
    if not files:
        print("未找到任何 results_*.csv 文件！")
        return
    
    print("正在扫描所有文件以获取 mdot 全局范围...")
    mdot_min = float('inf')
    mdot_max = -float('inf')
    valid_data = []
    for f in files:
        df = load_data(f)
        if df is None:
            continue
        valid_data.append((f, df))
        mdot_min = min(mdot_min, df['mdot'].min())
        mdot_max = max(mdot_max, df['mdot'].max())
    
    if not valid_data:
        print("没有有效的数据文件，退出。")
        return
    
    if mdot_max <= mdot_min:
        mdot_max = mdot_min + 1e-6
    print(f"全局 mdot 范围: {mdot_min:.4e} ~ {mdot_max:.4e}")
    
    print(f"找到 {len(valid_data)} 个有效文件，开始生成 mdot 动画帧...")
    frame_files = []
    for idx, (file_path, df) in enumerate(valid_data):
        step_num = re.findall(r'\d+', os.path.basename(file_path))[0]
        out_path = os.path.join(temp_dir, f'frame_{idx:04d}.png')
        try:
            draw_frame(df, step_num, out_path, vmin=mdot_min, vmax=mdot_max)
            frame_files.append(out_path)
        except Exception as e:
            print(f"绘制第 {step_num} 步时出错：{e}")
        if (idx+1) % 10 == 0:
            print(f"已生成 {idx+1}/{len(valid_data)} 帧")
    
    if not frame_files:
        print("没有成功生成任何帧，退出。")
        return
    
    print("所有帧生成完毕，开始合成 GIF...")
    images = [imageio.imread(f) for f in frame_files]
    gif_path = 'mdot_evolution.gif'
    imageio.mimsave(gif_path, images, fps=5, loop=0)
    print(f"mdot 动画已保存为 {gif_path}")

if __name__ == "__main__":
    main()