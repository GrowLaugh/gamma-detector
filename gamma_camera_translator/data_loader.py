import numpy as np
import os
import glob
import tensorflow as tf
from sklearn.model_selection import train_test_split

def load_data(data_dir, target_shape=(32, 32)):
    """
    加载数据，将输入(sim) resize 到与输出(real)相同的尺寸，并进行归一化。
    """
    sim_files = sorted(glob.glob(os.path.join(data_dir, 'sim', '*.npy')))
    real_files = sorted(glob.glob(os.path.join(data_dir, 'real', '*.npy')))
    
    if not sim_files or not real_files:
        raise FileNotFoundError(f"在 '{data_dir}' 中未找到 .npy 文件。请先运行 generate_dummy_data.py")
        
    if len(sim_files) != len(real_files):
        raise ValueError("sim 和 real 目录中的文件数量不匹配!")

    X = [] # Sim data (Input)
    Y = [] # Real data (Target)

    for sim_file, real_file in zip(sim_files, real_files):
        sim_data = np.load(sim_file)
        real_data = np.load(real_file)
        
        X.append(sim_data)
        Y.append(real_data)

    # 转换为Numpy数组
    X = np.array(X, dtype=np.float32)
    Y = np.array(Y, dtype=np.float32)

    # --- 关键预处理 ---
    
    # 1. 将输入(Sim) resize 到目标尺寸 (e.g., 60x60 -> 32x32)
    # 我们需要先增加一个通道维度
    X = X[..., tf.newaxis]
    X_resized = tf.image.resize(X, target_shape).numpy()
    
    # 2. 增加目标(Real)的通道维度
    Y = Y[..., tf.newaxis]

    # 3. 归一化 (使用训练集的最大值)
    # 注意：在真实项目中，应保存训练集的最大值，并用它来归一化验证集和测试集
    X_max = X_resized.max()
    Y_max = Y.max()
    
    if X_max == 0: X_max = 1.0
    if Y_max == 0: Y_max = 1.0
    
    X_norm = X_resized / X_max
    Y_norm = Y / Y_max
    
    print(f"数据加载完成。 X_shape: {X_norm.shape}, Y_shape: {Y_norm.shape}")
    
    # 4. 划分训练集和验证集
    X_train, X_val, Y_train, Y_val = train_test_split(
        X_norm, Y_norm, test_size=0.2, random_state=42
    )
    
    # 保存归一化参数以便预测时使用
    norm_params = {'X_max': X_max, 'Y_max': Y_max}
    
    return (X_train, Y_train), (X_val, Y_val), norm_params
