% =========================================================================
% === 脚本: 将 3D 参考库合成为 2D 参考库 (基于物理衰减加权) ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入: 3D 参考库文件 ---
input_3d_library_file = 'E:\110work\library生成\new48.5mmstep2.5mm.mat';

% --- 输出: 2D 参考库文件 ---
output_2d_library_file = 'E:\110work\library生成\projectionslibrary_symmetric_full2D2.5mm_Synthetic.mat';

% --- 物理参数 (NaI @ 140 keV) ---
mu_cm = 3.098;         % 线性衰减系数 (cm^-1)
mu_mm = mu_cm / 10;    % 转换为 mm^-1 (约 0.3098)

%% 2. 加载 3D 参考库
% -------------------------------------------------------------------------
fprintf('正在加载 3D 参考库: %s ...\n', input_3d_library_file);
if ~exist(input_3d_library_file, 'file')
    error('错误: 找不到输入文件。');
end

% 加载所需的变量
% 必须包含: lightMapLibrary (原始光分布图), x_coords, y_coords, z_coords
loaded_data = load(input_3d_library_file);

if ~isfield(loaded_data, 'lightMapLibrary')
    error('错误: 3D 库文件中缺少 lightMapLibrary 变量 (原始光分布矩阵)。无法进行合成。');
end

lightMap3D = loaded_data.lightMapLibrary;
x_coords = loaded_data.x_coords;
y_coords = loaded_data.y_coords;
z_coords = loaded_data.z_coords;

[nx, ny, nz] = size(lightMap3D);
fprintf('加载成功！库尺寸: %d(X) x %d(Y) x %d(Z)\n', nx, ny, nz);

%% 3. 计算深度权重 (Weights Calculation)
% -------------------------------------------------------------------------
% 假设 z_coords 代表物理深度 (或者距离表面的距离)
% 如果 z_coords 是负值 (例如 -0.5, -1.5)，取绝对值作为穿透深度
depths_mm = abs(z_coords);

% 1. 计算原始指数权重: W = exp(-mu * depth)
raw_weights = exp(-mu_mm * depths_mm);

% 2. 归一化权重 (使总和为 1)
norm_weights = raw_weights / sum(raw_weights);

% --- 打印权重信息供检查 ---
fprintf('============================================================\n');
fprintf('深度加权信息 (mu = %.4f mm^-1):\n', mu_mm);
fprintf('------------------------------------------------------------\n');
fprintf('  层索引 |  物理深度(mm) |  原始权重  |  归一化权重 (%%)\n');
for i = 1:nz
    fprintf('    %d    |    %6.2f     |   %.4f   |    %5.2f%%\n', ...
            i, depths_mm(i), raw_weights(i), norm_weights(i)*100);
end
fprintf('============================================================\n');

% 可视化权重分布
figure('Name', 'Depth Weighting');
bar(depths_mm, norm_weights);
xlabel('Depth (mm)'); ylabel('Normalized Weight');
title('Z-Axis Weighting based on Attenuation');
grid on;

%% 4. 执行合成 (3D -> 2D)
% -------------------------------------------------------------------------
fprintf('开始合成 2D 参考库...\n');

% 初始化新的 2D 容器
referprojection_2D = cell(nx, ny);   % 存储归一化投影向量 (用于定位)
lightMapLibrary_2D = cell(nx, ny);   % 存储合成后的光分布图 (用于显示)

% 获取探测器像素尺寸 (用于初始化空矩阵)
% 假设所有光分布图尺寸一致，取第一个非空的样本
first_valid = find(~cellfun(@isempty, lightMap3D), 1);
if isempty(first_valid), error('库是空的！'); end
sample_map = lightMap3D{first_valid};
map_size = size(sample_map);

% 循环遍历每个 (X, Y) 位置
count = 0;
total = nx * ny;

for ix = 1:nx
    for iy = 1:ny
        count = count + 1;
        
        % --- 核心合成逻辑 ---
        weighted_sum_map = zeros(map_size); % 初始化累加器
        has_data = false;
        
        for iz = 1:nz
            current_map = lightMap3D{ix, iy, iz};
            
            if ~isempty(current_map)
                % 转换为 double 进行计算
                current_map = double(current_map);
                
                % 加权累加: Final = Sum( Map_i * Weight_i )
                weighted_sum_map = weighted_sum_map + (current_map * norm_weights(iz));
                has_data = true;
            end
        end
        
        if has_data
            % 1. 存储合成后的光分布图
            lightMapLibrary_2D{ix, iy} = weighted_sum_map;
            
            % 2. 计算投影向量 (X 和 Y 方向)
            x_proj = sum(weighted_sum_map, 1); % 列求和 -> X投影
            y_proj = sum(weighted_sum_map, 2)'; % 行求和 -> Y投影 (转置为行向量)
            
            % 3. 生成归一化的参考投影 (用于 simi/roughest)
            combined_proj = [x_proj, y_proj];
            total_counts = sum(combined_proj);
            
            if total_counts > 0
                % 归一化到和为 1
                norm_proj = combined_proj / total_counts;
            else
                norm_proj = zeros(size(combined_proj));
            end
            
            referprojection_2D{ix, iy} = norm_proj;
        else
            % 如果该位置没有任何数据
            lightMapLibrary_2D{ix, iy} = [];
            referprojection_2D{ix, iy} = [];
        end
        
        if mod(count, 500) == 0
            fprintf('进度: %.1f%%\n', count/total*100);
        end
    end
end

%% 5. 保存结果
% -------------------------------------------------------------------------
fprintf('正在保存 2D 参考库: %s ...\n', output_2d_library_file);

% 重命名变量以匹配标准格式 (方便后续脚本直接 load)
referprojection = referprojection_2D;
lightMapLibrary = lightMapLibrary_2D;

% 保存 (只保留 x_coords 和 y_coords，z_coords 不再需要或设为 [])
% 也可以保留 z_coords 说明这是从哪些深度合成的
save(output_2d_library_file, 'referprojection', 'lightMapLibrary', 'x_coords', 'y_coords', 'z_coords', 'norm_weights');

fprintf('============================================================\n');
fprintf('2D 参考库生成完毕！\n');
fprintf('  - 原始尺寸: %d x %d x %d\n', nx, ny, nz);
fprintf('  - 合成尺寸: %d x %d\n', nx, ny);
fprintf('您可以直接在 step3 脚本中使用此文件进行 2D 定位。\n');