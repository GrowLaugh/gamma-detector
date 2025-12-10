% =========================================================================
% === 脚本 3 (批量版): 使用【混合ML + 3D插值算法】批量处理并汇总显示 ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入文件设置 ---
input_data_folder   = 'calib_data';       % 存放 calib_events*.mat 文件的文件夹
file_basename       = 'calib_events_p%04d.mat'; % 文件名格式

% --- !!! 在此处定义所有需要处理的文件编号 !!! ---
file_indices_to_process = [ ...
    503:508, ...  % 第一批
    603:608, ...  % 第二批
    703:708, ...  % 第三批
    803:808  ...  % 第四批
];

% --- 参考库设置 ---
% 请修改为您实际的3D参考库路径
libraryFile           = 'E:\110work\library生成\projectionslibrary_symmetric_full3D2.5mm.mat'; % 您的3D参考库

% --- 算法参数 ---
coarse_radius    = 5;  % 粗定位邻域半径
interp_padding   = 1;  % 插值边界
detector_size_mm = 50; % 探测器物理尺寸
display_res      = 1000; % 显示分辨率

%% 2. 加载参考库并初始化
% -------------------------------------------------------------------------
fprintf('正在加载参考库: %s ...\n', libraryFile);
if ~exist(libraryFile, 'file')
    error('错误: 找不到参考库文件 %s。请检查路径。', libraryFile);
end
load(libraryFile); % 加载 referprojection, x_coords, y_coords, z_coords
[libx, liby, libz] = size(referprojection);
fprintf('参考库加载成功！尺寸: %d x %d x %d\n', libx, liby, libz);

num_files = length(file_indices_to_process);
fprintf('准备批量处理 %d 个文件...\n', num_files);

% 用于存储每个点的最终均值位置 (X, Y, Z)
all_mean_positions = nan(num_files, 3); 

%% 3. 初始化绘图窗口和画布
% -------------------------------------------------------------------------
h_fig = figure('Name', '批量ML定位结果汇总', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);
axis_limit = detector_size_mm / 2;
display_edges = linspace(-axis_limit, axis_limit, display_res + 1);
display_centers = (display_edges(1:end-1) + display_edges(2:end)) / 2;
composite_image = zeros(display_res, display_res); % 主画布

% 启动并行池 (如果尚未启动)
if isempty(gcp('nocreate'))
    parpool; 
end

%% 4. 主循环：逐个处理文件
% -------------------------------------------------------------------------
fprintf('============================================================\n');
tic; % 开始计时

loop_counter = 0;
for file_idx = file_indices_to_process
    loop_counter = loop_counter + 1;
    
    % --- a. 加载当前文件 ---
    current_filename = sprintf(file_basename, file_idx);
    inputFile = fullfile(input_data_folder, current_filename);
    
    if ~exist(inputFile, 'file')
        fprintf('警告: 找不到文件 %s，已跳过。\n', current_filename);
        continue;
    end
    
    fprintf('------------------------------------------------------------\n');
    fprintf('正在处理文件 (%d/%d): %s\n', loop_counter, num_files, current_filename);
    
    data_struct = load(inputFile); 
    planeset = data_struct.planeset;
    pixels_xy = data_struct.pixels_xy;
    total_events = length(planeset);
    
    if total_events == 0
        fprintf('  警告: 文件中没有有效事件。\n');
        continue;
    end
    
    % --- b. 对当前文件执行 ML + 3D 插值算法 (并行计算) ---
    % 预分配结果矩阵 (N x 3)
    file_calc_indices = zeros(total_events, 3); 
    
    parfor i = 1:total_events
        % 1. 准备投影
        lightMap = planeset{i};
        x_proj = sum(lightMap, 1);
        y_proj = sum(lightMap, 2)';
        norm_proj = [x_proj, y_proj] / (sum(lightMap(:)) + eps);
        
        % 2. 粗定位
        [xrough, yrough] = roughest(norm_proj, pixels_xy, 3, 2, libx, liby);
        
        % 3. 确定计算邻域
        [~, rx_s, ry_s] = neighbourdeter(xrough, yrough, referprojection, coarse_radius);
        
        rx_calc_start = max(1, rx_s(1) - interp_padding);
        rx_calc_end   = min(libx, rx_s(end) + interp_padding);
        rx_calc = rx_calc_start:rx_calc_end;
        
        ry_calc_start = max(1, ry_s(1) - interp_padding);
        ry_calc_end   = min(liby, ry_s(end) + interp_padding);
        ry_calc = ry_calc_start:ry_calc_end;
        
        % 4. 3D ML 搜索
        sz_x = length(rx_calc); sz_y = length(ry_calc);
        L_sub = ones(sz_x, sz_y, libz) * -inf;
        
        % 注意：在 parfor 中不能直接访问外部切片 cell，最好传递整个 cell 或特定部分
        % 这里为了简化，直接访问 referprojection (Matlab 会自动处理广播变量)
        for xi = 1:sz_x
            for yi = 1:sz_y
                xg = rx_calc(xi); yg = ry_calc(yi);
                for z = 1:libz
                    ref = referprojection{xg, yg, z};
                    if ~isempty(ref)
                        L_sub(xi, yi, z) = calc_logL_gaussian(norm_proj, ref);
                    end
                end
            end
        end
        
        % 5. 寻找峰值并插值
        % 为了避免边缘效应，在原始搜索范围内找最大值
        x_s_start = rx_s(1) - rx_calc_start + 1;
        x_s_end   = rx_s(end) - rx_calc_start + 1;
        y_s_start = ry_s(1) - ry_calc_start + 1;
        y_s_end   = ry_s(end) - ry_calc_start + 1;
        
        L_search = ones(size(L_sub)) * -inf;
        L_search(x_s_start:x_s_end, y_s_start:y_s_end, :) = ...
            L_sub(x_s_start:x_s_end, y_s_start:y_s_end, :);
            
        [max_L, lin_idx] = max(L_search(:));
        
        if max_L > -inf
            [sx, sy, sz] = ind2sub(size(L_search), lin_idx);
            [ix, iy, iz] = interpolate_peak_parabolic(L_sub, sx, sy, sz);
            
            % 转换回全局浮点索引
            file_calc_indices(i, :) = [ix + rx_calc_start - 1, ...
                                       iy + ry_calc_start - 1, ...
                                       iz];
        else
            file_calc_indices(i, :) = [NaN, NaN, NaN];
        end
    end
    
    % --- c. 转换为物理坐标 ---
    valid_mask = ~isnan(file_calc_indices(:,1));
    valid_indices = file_calc_indices(valid_mask, :);
    
    if isempty(valid_indices)
        fprintf('  无有效定位结果。\n');
        continue;
    end
    
    % 索引 -> 物理坐标 (假设坐标是均匀分布的)
    % x_coords, y_coords, z_coords 来自加载的库
    dx = x_coords(2)-x_coords(1); x0 = x_coords(1);
    dy = y_coords(2)-y_coords(1); y0 = y_coords(1);
    % dz = z_coords(2)-z_coords(1); z0 = z_coords(1);
    
    pos_x_mm = (valid_indices(:,1) - 1) * dx + x0;
    pos_y_mm = (valid_indices(:,2) - 1) * dy + y0;
    % pos_z_mm = (valid_indices(:,3) - 1) * dz + z0;
    
    % --- d. 分析与累加 ---
    mean_pos = mean([pos_x_mm, pos_y_mm], 1);
    all_mean_positions(loop_counter, 1:2) = mean_pos;
    
    fprintf('  -> 均值位置: (%.2f, %.2f) mm, 有效点数: %d\n', ...
        mean_pos(1), mean_pos(2), size(valid_indices, 1));
    
    % 生成当前文件的 PSF 并累加到主画布
    current_psf = histcounts2(pos_y_mm, pos_x_mm, display_edges, display_edges);
    if max(current_psf(:)) > 0
        current_psf = current_psf / max(current_psf(:)); % 归一化，使每个点源亮度一致
    end
    composite_image = composite_image + current_psf;
    
    % --- e. 实时绘图 ---
    figure(h_fig);
    imagesc(display_centers, display_centers, composite_image);
    set(gca, 'YDir', 'normal');
    colormap('jet'); colorbar; axis equal;
    xlim([-axis_limit, axis_limit]); ylim([-axis_limit, axis_limit]);
    xlabel('X (mm)'); ylabel('Y (mm)');
    title(sprintf('批量处理中... (%d/%d)', loop_counter, num_files));
    
    % 添加白色网格 (5mm 间距)
    hold on;
    grid_ticks = -axis_limit:5:axis_limit;
    set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks);
    grid on;
    set(gca, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
    
    % 标记编号
    for k = 1:loop_counter
        if ~isnan(all_mean_positions(k, 1))
            text(all_mean_positions(k,1), all_mean_positions(k,2), ...
                num2str(file_indices_to_process(k)), ...
                'Color', 'w', 'FontSize', 8, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end
    hold off;
    drawnow;
end

toc;
fprintf('============================================================\n');
fprintf('处理完成！\n');

%% 本地函数
% =========================================================================
function [x_interp, y_interp, z_interp] = interpolate_peak_parabolic(L_matrix, x_idx, y_idx, z_idx)
% 3D 抛物线插值 (来自您的 MLE 脚本)
    [sub_libx, sub_liby, sub_libz] = size(L_matrix);
    
    % X
    xm1 = -inf; xp1 = -inf; xc = L_matrix(x_idx, y_idx, z_idx);
    if x_idx > 1, xm1 = L_matrix(x_idx-1, y_idx, z_idx); end
    if x_idx < sub_libx, xp1 = L_matrix(x_idx+1, y_idx, z_idx); end
    denom_x = 2 * (xm1 - 2*xc + xp1);
    if denom_x < 0 && isfinite(xm1) && isfinite(xp1)
        x_interp = x_idx + (xm1 - xp1) / denom_x;
    else
        x_interp = x_idx;
    end
    
    % Y
    ym1 = -inf; yp1 = -inf; yc = L_matrix(x_idx, y_idx, z_idx);
    if y_idx > 1, ym1 = L_matrix(x_idx, y_idx-1, z_idx); end
    if y_idx < sub_liby, yp1 = L_matrix(x_idx, y_idx+1, z_idx); end
    denom_y = 2 * (ym1 - 2*yc + yp1);
    if denom_y < 0 && isfinite(ym1) && isfinite(yp1)
        y_interp = y_idx + (ym1 - yp1) / denom_y;
    else
        y_interp = y_idx;
    end
    
    % Z
    zm1 = -inf; zp1 = -inf; zc = L_matrix(x_idx, y_idx, z_idx);
    if z_idx > 1, zm1 = L_matrix(x_idx, y_idx, z_idx-1); end
    if z_idx < sub_libz, zp1 = L_matrix(x_idx, y_idx, z_idx+1); end
    denom_z = 2 * (zm1 - 2*zc + zp1);
    if denom_z < 0 && isfinite(zm1) && isfinite(zp1)
        z_interp = z_idx + (zm1 - zp1) / denom_z;
    else
        z_interp = z_idx;
    end
end

function L = calc_logL_gaussian(norm_sample_proj, norm_ref_proj)
% ML 核心函数 (高斯近似)
    pseudolightyield = 8000; 
    norm_ref_proj_safe = norm_ref_proj + eps;
    sample_photons = norm_sample_proj * pseudolightyield;
    ref_photons_safe = norm_ref_proj_safe * pseudolightyield;
    
    term1 = (sample_photons - ref_photons_safe).^2 ./ (2 * ref_photons_safe);
    term2 = 0.5 * log(2 * pi * ref_photons_safe);
    
    L = sum(-term1 - term2);
end