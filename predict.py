import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
import sys
import os
import json

def predict_new_sample(model_path, params_path, sim_file_path, output_path='prediction.png'):
    """
    加载模型和新数据，进行预测并保存结果。
    """
    # --- 1. 加载模型和归一化参数 ---
    print(f"正在加载模型: {model_path}")
    model = tf.keras.models.load_model(model_path)
    
    print(f"正在加载归一化参数: {params_path}")
    with open(params_path, 'r') as f:
        norm_params = json.load(f)
    X_max = norm_params['X_max']
    Y_max = norm_params['Y_max']

    # --- 2. 加载和预处理新样本 ---
    print(f"正在加载输入样本: {sim_file_path}")
    sim_data = np.load(sim_file_path).astype(np.float32)
    
    # 应用与训练时完全相同的预处理
    # (a) 增加批次和通道维度
    input_data = sim_data[np.newaxis, ..., np.newaxis]
    # (b) Resize
    input_data_resized = tf.image.resize(input_data, (32, 32))
    # (c) 归一化
    input_data_norm = input_data_resized / X_max
    
    print("预处理完成，开始预测...")
    # --- 3. 预测 ---
    predicted_norm = model.predict(input_data_norm)
    
    # --- 4. 后处理 (反归一化) ---
    predicted_real_scale = predicted_norm * Y_max
    
    # 移除批次和通道维度
    predicted_image = predicted_real_scale.squeeze()
    
    print("预测完成。")
    # --- 5. 可视化和保存 ---
    plt.figure(figsize=(12, 6))
    
    # (a) 原始 60x60 输入
    plt.subplot(1, 2, 1)
    plt.imshow(sim_data, cmap='viridis')
    plt.title("原始输入 (Sim 60x60)")
    plt.axis('off')
    
    # (b) 预测的 32x32 输出
    plt.subplot(1, 2, 2)
    plt.imshow(predicted_image, cmap='inferno')
    plt.title("预测输出 (Real 32x32)")
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path)
    print(f"预测结果图已保存到: {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python predict.py <sim_data_file.npy>")
        print("例如: python predict.py data/sim/pos_001.npy")
        sys.exit(1)
        
    sim_file = sys.argv[1]
    if not os.path.exists(sim_file):
        print(f"错误: 文件 {sim_file} 不存在。")
        sys.exit(1)

    predict_new_sample(
        model_path='sim_to_real_model.h5',
        params_path='norm_params.json',
        sim_file_path=sim_file,
        output_path='new_prediction.png'
    )
