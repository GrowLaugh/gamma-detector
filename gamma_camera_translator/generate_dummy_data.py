import numpy as np
import os
from scipy.stats import multivariate_normal

def create_gaussian_blob(size, x0, y0, sigma=5):
    """创建一个高斯斑点"""
    x = np.arange(size[0])
    y = np.arange(size[1])
    X, Y = np.meshgrid(x, y)
    pos = np.empty(X.shape + (2,))
    pos[:, :, 0] = X
    pos[:, :, 1] = Y
    rv = multivariate_normal([x0, y0], [[sigma, 0], [0, sigma]])
    return rv.pdf(pos)

def generate_data(num_samples=50, data_path='data'):
    """生成模拟的 sim (60x60) 和 real (32x32) 数据"""
    sim_path = os.path.join(data_path, 'sim')
    real_path = os.path.join(data_path, 'real')
    os.makedirs(sim_path, exist_ok=True)
    os.makedirs(real_path, exist_ok=True)

    print(f"正在生成 {num_samples} 组虚拟数据...")
    for i in range(num_samples):
        # 随机中心点
        x_sim, y_sim = np.random.randint(10, 50, 2)
        x_real, y_real = (x_sim / 60) * 32, (y_sim / 60) * 32

        # 1. 创建 60x60 的“模拟”数据 (更清晰)
        sim_data = create_gaussian_blob((60, 60), x_sim, y_sim, sigma=5)
        sim_data = sim_data / sim_data.max() # 归一化
        
        # 2. 创建 32x32 的“真实”数据 (更模糊，有噪声)
        real_data = create_gaussian_blob((32, 32), x_real, y_real, sigma=7)
        real_data = real_data / real_data.max()
        real_data += np.random.rand(32, 32) * 0.1 # 添加噪声
        
        # 保存为 .npy 文件
        np.save(os.path.join(sim_path, f'pos_{i+1:03d}.npy'), sim_data.astype(np.float32))
        np.save(os.path.join(real_path, f'pos_{i+1:03d}.npy'), real_data.astype(np.float32))

    print(f"数据生成完毕，已保存至 '{data_path}' 目录。")

if __name__ == "__main__":
    generate_data(num_samples=200) # 生成200个训练样本
