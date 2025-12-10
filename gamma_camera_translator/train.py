import tensorflow as tf
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
import matplotlib.pyplot as plt
import numpy as np
import json

from data_loader import load_data
from model import build_unet

# --- 1. 设置参数 ---
DATA_DIR = 'data'
MODEL_SAVE_PATH = 'sim_to_real_model.h5'
PARAMS_SAVE_PATH = 'norm_params.json'
INPUT_SHAPE = (32, 32, 1)
EPOCHS = 100
BATCH_SIZE = 16

# --- 2. 加载数据 ---
print("正在加载和预处理数据...")
(X_train, Y_train), (X_val, Y_val), norm_params = load_data(DATA_DIR, target_shape=INPUT_SHAPE[:2])

# 保存归一化参数
with open(PARAMS_SAVE_PATH, 'w') as f:
    json.dump({k: float(v) for k, v in norm_params.items()}, f)
print(f"归一化参数已保存到 {PARAMS_SAVE_PATH}")

# --- 3. 构建和编译模型 ---
print("正在构建模型...")
model = build_unet(input_shape=INPUT_SHAPE)
model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
              loss='mean_squared_error', # MSE 是图像回归的标准损失
              metrics=['mean_absolute_error'])
model.summary()

# --- 4. 设置回调函数 ---
callbacks = [
    EarlyStopping(patience=10, verbose=1, restore_best_weights=True),
    ModelCheckpoint(MODEL_SAVE_PATH, verbose=1, save_best_only=True)
]

# --- 5. 训练模型 ---
print("开始训练模型...")
history = model.fit(
    X_train, Y_train,
    validation_data=(X_val, Y_val),
    epochs=EPOCHS,
    batch_size=BATCH_SIZE,
    callbacks=callbacks
)

print(f"训练完成。最佳模型已保存到 {MODEL_SAVE_PATH}")

# --- 6. 绘制训练历史 ---
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(history.history['loss'], label='训练集损失 (Loss)')
plt.plot(history.history['val_loss'], label='验证集损失 (Val Loss)')
plt.title('模型损失')
plt.xlabel('Epochs')
plt.ylabel('Loss (MSE)')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(history.history['mean_absolute_error'], label='训练集 MAE')
plt.plot(history.history['val_mean_absolute_error'], label='验证集 MAE')
plt.title('平均绝对误差 (MAE)')
plt.xlabel('Epochs')
plt.ylabel('MAE')
plt.legend()
plt.tight_layout()
plt.savefig('training_history.png')
print("训练历史图已保存到 training_history.png")

# --- 7. 可视化一些验证结果 ---
print("正在生成验证结果对比图...")
model.load_weights(MODEL_SAVE_PATH) # 确保加载最佳模型
Y_pred = model.predict(X_val)

n_samples = 5 # 显示5个样本
plt.figure(figsize=(15, 5 * n_samples))
for i in range(n_samples):
    # (a) 输入的模拟数据 (resize后)
    plt.subplot(n_samples, 3, i * 3 + 1)
    plt.imshow(X_val[i].squeeze(), cmap='viridis')
    plt.title(f"输入 (Sim) {i+1}")
    plt.axis('off')

    # (b) 预测的校准数据
    plt.subplot(n_samples, 3, i * 3 + 2)
    plt.imshow(Y_pred[i].squeeze(), cmap='inferno')
    plt.title(f"预测 (Predicted Real) {i+1}")
    plt.axis('off')
    
    # (c) 真实的校准数据
    plt.subplot(n_samples, 3, i * 3 + 3)
    plt.imshow(Y_val[i].squeeze(), cmap='inferno')
    plt.title(f"真实 (Ground Truth) {i+1}")
    plt.axis('off')

plt.tight_layout()
plt.savefig('validation_results.png')
print("验证结果对比图已保存到 validation_results.png")

