import cv2
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from scipy.spatial import cKDTree

def generate_jet_reference():
    """
    生成 matplotlib 的 jet(256) 颜色查找表。
    返回: 
        colormap_lut: (256, 3) 的数组，包含 jet 从 0 到 255 对应的 RGB 值
    """
    # 获取 matplotlib 的 jet 色标对象
    cmap = plt.get_cmap('jet')
    # 生成 0 到 255 的索引
    indices = np.linspace(0, 1, 256)
    # 获取对应的 RGB 值 (去掉了 Alpha 通道)，并扩展到 0-255 范围
    # matplotlib 返回的是 0-1 的 float，需要乘 255 转为 uint8
    colormap_lut = (cmap(indices)[:, :3] * 255).astype(np.uint8)
    return colormap_lut

def image_to_matrix(img_bgr, tree):
    """
    将图片反演回数值矩阵。
    img_bgr: 输入的 BGR 图片
    tree: 用于快速查找颜色的 KDTree
    """
    # 1. OpenCV 读入是 BGR，Matplotlib 的 jet 是 RGB，需要转换
    img_rgb = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB)
    
    # 2. 获取图片尺寸
    h, w, _ = img_rgb.shape
    
    # 3. 将图片展平为 (N, 3) 的数组，以便进行批量查找
    pixels = img_rgb.reshape(-1, 3)
    
    # 4. 核心步骤：查找每个像素在 jet 色标中最近的索引值 (0-255)
    # dist 是距离，idx 是查找到的索引 (即我们要的相对能量值)
    # k=1 表示只找最近的那一个
    _, idx = tree.query(pixels, k=1)
    
    # 5. 将一维索引还原回 (H, W) 的矩阵
    value_matrix = idx.reshape(h, w)
    
    return value_matrix

def main(input_folder, output_filename="final_accumulated_heatmap.png"):
    # 1. 准备 Jet 色标查找树
    print("正在生成 Jet 色标查找表...")
    jet_lut = generate_jet_reference()
    # 使用 KDTree 建立快速索引，加速颜色匹配过程
    tree = cKDTree(jet_lut)

    # 2. 获取图片
    image_paths = []
    for ext in ['*.jpg', '*.jpeg', '*.png', '*.bmp', '*.tif']:
        image_paths.extend(glob.glob(os.path.join(input_folder, ext)))
    image_paths.sort()

    if not image_paths:
        print("未找到图片。")
        return

    print(f"找到 {len(image_paths)} 张图片。目标尺寸: 1000x1000")

    # 3. 初始化累加矩阵
    # 使用 float64 防止累加溢出
    target_size = (1000, 1000)
    accumulated_matrix = np.zeros(target_size, dtype=np.float64)

    for i, img_path in enumerate(image_paths):
        img = cv2.imread(img_path)
        if img is None: continue

        # (1) 调整大小到 1000x1000
        # 使用 INTER_CUBIC 或 LANCZOS4 保证缩放平滑
        resized_img = cv2.resize(img, target_size, interpolation=cv2.INTER_CUBIC)

        # (2) 反演回数值矩阵 (Matrix Recovery)
        # 这一步会将 RGB 颜色转回 0-255 的能量值
        current_matrix = image_to_matrix(resized_img, tree)
        
        # (3) 累加
        accumulated_matrix += current_matrix
        
        print(f"[{i+1}/{len(image_paths)}] 已叠加: {os.path.basename(img_path)}")

    # 4. 生成最终结果
    if len(image_paths) > 0:
        print("正在生成最终热力图...")
        
        # 创建画布，不带边框
        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = plt.axes([0, 0, 1, 1]) # 占满整个画布
        ax.set_axis_off()
        
        # 绘制热力图
        # cmap='jet': 再次使用 jet 上色
        # accumulation 矩阵中的大数值会自动对应红色，小数值(背景)对应蓝色
        plt.imshow(accumulated_matrix, cmap='jet', interpolation='nearest')
        
        # 保存
        plt.savefig(output_filename, bbox_inches='tight', pad_inches=0, dpi=100)
        plt.close()
        
        print("-" * 30)
        print(f"处理完成！原理：RGB -> 数值 -> 累加 -> RGB")
        print(f"结果已保存: {output_filename}")
        
        # (可选) 如果你需要保存累加后的原始数据以供后续分析
        # np.save("accumulated_data.npy", accumulated_matrix)
        # print("原始数据矩阵已保存为 accumulated_data.npy")

if __name__ == "__main__":
    input_dir = "./gamma_images" # 图片文件夹路径
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)
        print(f"请将图片放入 '{input_dir}'")
    else:
        main(input_dir)