import cv2
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from scipy.spatial import cKDTree

def generate_jet_reference():
    """
    生成 matplotlib 的 jet(256) 颜色查找表。
    """
    cmap = plt.get_cmap('jet')
    indices = np.linspace(0, 1, 256)
    colormap_lut = (cmap(indices)[:, :3] * 255).astype(np.uint8)
    return colormap_lut

def image_to_matrix(img_bgr, tree):
    """
    将图片反演回数值矩阵。
    """
    img_rgb = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB)
    h, w, _ = img_rgb.shape
    pixels = img_rgb.reshape(-1, 3)
    _, idx = tree.query(pixels, k=1)
    value_matrix = idx.reshape(h, w)
    return value_matrix

def main(input_folder, output_filename="final_cropped_heatmap.png"):
    # 1. 准备 Jet 色标查找树
    print("正在生成 Jet 色标查找表...")
    jet_lut = generate_jet_reference()
    tree = cKDTree(jet_lut)

    # 2. 获取图片
    image_paths = []
    for ext in ['*.jpg', '*.jpeg', '*.png', '*.bmp', '*.tif']:
        image_paths.extend(glob.glob(os.path.join(input_folder, ext)))
    image_paths.sort()

    if not image_paths:
        print("未找到图片。")
        return

    # --- 设置参数 ---
    resize_target = (1000, 1000) # 先缩放到 1000x1000
    
    # 定义裁剪比例 (10% 到 90%)
    # 对应坐标: 100 到 900
    crop_start_ratio = 0.1 
    crop_end_ratio = 0.9   
    
    print(f"找到 {len(image_paths)} 张图片。")
    print(f"处理逻辑: 缩放至 {resize_target} -> 裁剪区间 {crop_start_ratio*100}%-{crop_end_ratio*100}%")

    # 3. 初始化累加矩阵
    # 注意：累加矩阵的大小不再是 1000x1000，而是裁剪后的大小
    # 我们将在处理第一张图时动态初始化它
    accumulated_matrix = None

    for i, img_path in enumerate(image_paths):
        img = cv2.imread(img_path)
        if img is None: continue

        # (1) 调整大小到 1000x1000
        resized_img = cv2.resize(img, resize_target, interpolation=cv2.INTER_CUBIC)

        # (2) 计算裁剪坐标并执行裁剪 (新增步骤)
        h, w = resized_img.shape[:2]
        y_start = int(h * crop_start_ratio) # 100
        y_end = int(h * crop_end_ratio)     # 900
        x_start = int(w * crop_start_ratio) # 100
        x_end = int(w * crop_end_ratio)     # 900

        # Python 切片语法: [y_start:y_end, x_start:x_end]
        cropped_img = resized_img[y_start:y_end, x_start:x_end]

        # (3) 初始化累加矩阵 (仅在第一次运行时执行)
        if accumulated_matrix is None:
            crop_h, crop_w = cropped_img.shape[:2]
            accumulated_matrix = np.zeros((crop_h, crop_w), dtype=np.float64)
            print(f"裁剪后实际尺寸为: {crop_w}x{crop_h} (请确认是否符合预期)")

        # (4) 反演回数值矩阵 (针对裁剪后的图)
        current_matrix = image_to_matrix(cropped_img, tree)
        
        # (5) 累加
        accumulated_matrix += current_matrix
        
        print(f"[{i+1}/{len(image_paths)}] 已处理: {os.path.basename(img_path)}")

    # 4. 生成最终结果
    if accumulated_matrix is not None:
        print("正在生成最终热力图...")
        
        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = plt.axes([0, 0, 1, 1])
        ax.set_axis_off()
        
        plt.imshow(accumulated_matrix, cmap='jet', interpolation='nearest')
        
        plt.savefig(output_filename, bbox_inches='tight', pad_inches=0, dpi=100)
        plt.close()
        
        print("-" * 30)
        print(f"处理完成！")
        print(f"结果已保存: {output_filename}")

if __name__ == "__main__":
    input_dir = "./gamma_images" # 图片文件夹路径
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)
        print(f"请将图片放入 '{input_dir}'")
    else:
        main(input_dir)
