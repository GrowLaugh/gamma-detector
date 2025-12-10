% =========================================================================
% === 脚本 3 (批量版): 使用【2D MLE + 插值】批量处理 2D参考库数据 ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入文件设置 ---
input_data_folder   = 'E:\110work\ExpDataProc\calib_data1';
file_basename       = 'calib_events_p%04d.mat'; % 文件名格式 

% --- !!! 在此处定义所有需要处理的文件编号 !!! ---
file_indices_to_process = [ ...
    % 503:508, ... 
     603:608, ... 
    % 703:708, ... 
    % 803:808, ...  
    % 901, 902, 909 ...
];

% --- 输出设置 ---
output_image_folder = 'mle_Results1mm_calib1-2Dtest';   % 结果保存文件夹
if ~exist(output_image_folder, 'dir'), mkdir(output_image_folder); end

% --- 参考库设置 ---
% 您的 2D 参考库路径
libraryFile = 'E:\110work\library生成\projectionslibrary_symmetric_full2Dleijia25.mat';

% --- 物理坐标参数 ---
lib_phys_start        = -25.0; % 参考库起始坐标 (mm)
lib_phys_end          = 25.0;  % 参考库结束坐标 (mm)
detector_size_mm      = 50;    % 探测器物理尺寸 (mm)
pixel_size_mm         = 6.0625; % 探测器单像素尺寸 (48.5mm / 8 = 6.0625)

% --- 算法参数 ---
coarse_radius    = 5;  % 粗定位邻域半径
interp_padding   = 1;  % 插值边界

% --- 显示与网格参数 ---
display_resolution    = 512;    % 最终图像的分辨率 (像素)
grid_interval_mm      = 5;       % 白色网格线间距

% --- 似然矩阵平滑参数 ---
enable_ml_smoothing   = false;    % 设为 true 开启似然矩阵平滑 (抗网格化)
ml_smoothing_sigma    = 0.5;     % 平滑强度

% --- 结果图显示平滑 ---
enable_result_smoothing = false;  
result_smoothing_sigma  = 1.0;   

%% 2. 加载参考库 (只加载一次)
% -------------------------------------------------------------------------
fprintf('正在加载参考库: %s ...\n', libraryFile);
if ~exist(libraryFile, 'file')
    error('错误: 找不到参考库文件 %s', libraryFile);
end
load(libraryFile); % 加载 referprojection, x_coords, y_coords (2D库通常没有 z_coords)

% 检查参考库维度
[libx, liby, libz] = size(referprojection);
if libz > 1
    fprintf('警告: 检测到参考库有 Z 维度 (%d)，但本脚本将只使用第 1 层作为 2D 库。\n', libz);
    libz = 1; % 强制设为1，适配下面的逻辑
    % 如果需要，可以提取中间层: referprojection = referprojection(:,:,ceil(end/2));
else
    fprintf('参考库加载成功！尺寸: %d x %d (2D)\n', libx, liby);
end

% --- 初始化汇总画布 ---
axis_limit = detector_size_mm / 2; 
bin_edges = linspace(-axis_limit, axis_limit, display_resolution + 1);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
composite_image = zeros(display_resolution, display_resolution); % 主画布

num_files = length(file_indices_to_process);
fprintf('准备批量处理 %d 个文件...\n', num_files);

% 启动并行池
if isempty(gcp('nocreate')), parpool; end

%% 3. 主循环：逐个处理文件
% -------------------------------------------------------------------------
for k = 1:num_files
    file_idx = file_indices_to_process(k);
    
    % --- a. 加载当前文件 ---
    current_filename = sprintf(file_basename, file_idx);
    inputFile = fullfile(input_data_folder, current_filename);
    
    if ~exist(inputFile, 'file')
        fprintf('警告: 找不到文件 %s，已跳过。\n', current_filename);
        continue;
    end
    
    fprintf('------------------------------------------------------------\n');
    fprintf('正在处理 (%d/%d): %s ... ', k, num_files, current_filename);
    
    data_struct = load(inputFile);
    planeset = data_struct.planeset;
    pixels_xy = data_struct.pixels_xy;
    total_events = length(planeset);
    
    if total_events == 0
        fprintf('无有效事件。\n'); continue;
    end
    
    % --- b. 执行 2D ML + 插值定位算法 (并行) ---
    file_calc_indices = zeros(total_events, 2); % 只需要 X, Y
    
    parfor i = 1:total_events
        lightMap = planeset{i};
        x_proj = sum(lightMap, 1);
        y_proj = sum(lightMap, 2)';
        
        combined_proj = [x_proj, y_proj];
        norm_proj = combined_proj / (sum(combined_proj) + eps);
        
        % 1. 粗定位 (使用物理参数)
        try
            [xrough, yrough] = roughestnew(norm_proj, pixels_xy, 8, 2, ...
                                           libx, liby, pixel_size_mm, ...
                                           lib_phys_start, lib_phys_end);
        catch
            xrough = round(libx/2); yrough = round(liby/2);
        end
        
        % 2. 确定邻域 (2D)
        [~, rx_s, ry_s] = neighbourdeter_2d(xrough, yrough, referprojection, coarse_radius);
        
        rx_calc_start = max(1, rx_s(1) - interp_padding);
        rx_calc_end   = min(libx, rx_s(end) + interp_padding);
        rx_calc = rx_calc_start:rx_calc_end;
        
        ry_calc_start = max(1, ry_s(1) - interp_padding);
        ry_calc_end   = min(liby, ry_s(end) + interp_padding);
        ry_calc = ry_calc_start:ry_calc_end;
        
        % 3. 2D ML 搜索
        sz_x = length(rx_calc); 
        sz_y = length(ry_calc);
        L_sub = ones(sz_x, sz_y) * -inf; % 2D 似然矩阵

        
        for xi = 1:sz_x
            for yi = 1:sz_y
                xg = rx_calc(xi); yg = ry_calc(yi);
                % 获取参考投影 (注意处理 cell 维度)
                if ndims(referprojection) == 2
                    ref = referprojection{xg, yg};
                else
                    ref = referprojection{xg, yg, 1}; % 如果库是3D的，取第1层
                end
                
                if ~isempty(ref)
                    L_sub(xi, yi) = calc_logL_gaussian(norm_proj, ref);
                end
            end
        end
        
        % 4. 似然矩阵平滑 (2D)
        if enable_ml_smoothing
            valid_vals = L_sub(isfinite(L_sub));
            if ~isempty(valid_vals)
                min_val = min(valid_vals);
                L_sub(L_sub == -inf) = min_val - 10; 
                try
                    L_sub = imgaussfilt(L_sub, ml_smoothing_sigma);
                catch
                end
            end
        end
        
        % 5. 寻找峰值
        x_s_start = rx_s(1) - rx_calc_start + 1;
        x_s_end   = rx_s(end) - rx_calc_start + 1;
        y_s_start = ry_s(1) - ry_calc_start + 1;
        y_s_end   = ry_s(end) - ry_calc_start + 1;
        
        L_search = ones(size(L_sub)) * -inf;
        L_search(x_s_start:x_s_end, y_s_start:y_s_end) = ...
            L_sub(x_s_start:x_s_end, y_s_start:y_s_end);
            
        [max_L, lin_idx] = max(L_search(:));
        
        if max_L > -inf
            [sx, sy] = ind2sub(size(L_search), lin_idx);
            % 6. 2D 抛物线插值
            [ix, iy] = interpolate_peak_parabolic_2d(L_sub, sx, sy);
            
            % 转换回全局浮点索引
            file_calc_indices(i, :) = [ix + rx_calc_start - 1, ...
                                       iy + ry_calc_start - 1];
        else
            file_calc_indices(i, :) = [NaN, NaN];
        end
    end
    
    % 移除无效点
    valid_mask = ~isnan(file_calc_indices(:,1));
    file_calc_indices = file_calc_indices(valid_mask, :);
    
    if isempty(file_calc_indices)
        fprintf('定位失败 (无有效解)。\n'); continue;
    end
    
    % --- c. 转换为物理坐标 (mm) ---
    dx = x_coords(2)-x_coords(1); x0 = x_coords(1);
    dy = y_coords(2)-y_coords(1); y0 = y_coords(1);
    
    xpo_mm = (file_calc_indices(:,1) - 1) * dx + x0;
    ypo_mm = (file_calc_indices(:,2) - 1) * dy + y0;
    
    % --- d. 生成直方图数据 (修正：使用 histcounts2 替代 histogram2) ---
    % histcounts2 仅计算数据不绘图，速度快且不弹窗
    % 注意：输出矩阵的行对应Y (edge1)，列对应X (edge2)
    raw_counts = histcounts2(ypo_mm, xpo_mm, bin_edges, bin_edges);
    
    % 旋转/转置以匹配 imagesc 的坐标系 (X轴水平，Y轴垂直)
    % 通常 imagesc(x, y, C) 中 C(i,j) 对应 y(i) 和 x(j)
    % histcounts2 返回 N(i,j) 对应 edges1(i) 和 edges2(j) -> 这里是 Y 和 X
    % 所以 raw_counts 不需要转置即可直接给 imagesc(x, y, raw_counts) 使用?
    % 让我们保持之前的逻辑：
    current_img = raw_counts; 
    
    if enable_result_smoothing
        current_img = imgaussfilt(current_img, result_smoothing_sigma);
    end
    
    % ================= 保存图片部分 =================
    
    % --- 1. 保存 2D 结果图 (带白色网格) ---
    fig_2d = figure('Name', '2D Result', 'Visible', 'off');
    imagesc(bin_centers, bin_centers, current_img);
    set(gca, 'YDir', 'normal');
    colormap('jet'); colorbar; axis equal tight;
    title(sprintf('2D MLE: p%04d (N=%d)', file_idx, length(xpo_mm)));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % 添加白色网格
    hold on;
    grid_ticks = -axis_limit:grid_interval_mm:axis_limit;
    set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks);
    grid on;
    set(gca, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
    set(gca, 'Layer', 'top');
    hold off;
    
    exportgraphics(fig_2d, fullfile(output_image_folder, sprintf('MLE_2D_p%04d.png', file_idx)));
    close(fig_2d);
    
    % --- 2. 保存 3D 结果图 (使用 histogram2 绘图) ---
    fig_3d = figure('Name', '3D Result', 'Visible', 'off');
    histogram2(xpo_mm, ypo_mm, bin_edges, bin_edges, 'FaceColor', 'flat', 'ShowEmptyBins', 'off');
    view(45, 30); 
    colormap('jet'); colorbar;
    title(sprintf('3D MLE: p%04d', file_idx));
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Counts');
    
    exportgraphics(fig_3d, fullfile(output_image_folder, sprintf('MLE_3D_p%04d.png', file_idx)));
    close(fig_3d);
    
    % --- e. 累加到汇总图 ---
    if max(current_img(:)) > 0
        % 归一化后累加，使每个点源亮度一致
        composite_image = composite_image + (current_img / max(current_img(:)));
    end
    
    fprintf('已完成保存。\n');
end

fprintf('============================================================\n');

%% 4. 生成并保存最终汇总图 (带白色网格)
% -------------------------------------------------------------------------
fprintf('正在生成汇总图...\n');
h_sum = figure('Name', 'Summary', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);

imagesc(bin_centers, bin_centers, composite_image);
set(gca, 'YDir', 'normal');
colormap('jet'); colorbar; axis equal tight;

% --- 添加白色网格 ---
hold on;
grid_ticks = -axis_limit:grid_interval_mm:axis_limit; 
set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks);
grid on;
set(gca, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
set(gca, 'Layer', 'top'); 
hold off;

title(sprintf('MLE定位结果汇总 (%d 个点源)', num_files));
xlabel('X Position (mm)');
ylabel('Y Position (mm)');

summary_file = fullfile(output_image_folder, 'MLE_Summary_Combined.png');
exportgraphics(h_sum, summary_file, 'Resolution', 300);

fprintf('汇总图已保存至: %s\n', summary_file);
fprintf('任务全部完成！\n');

%% 本地函数 (MLE + 插值核心)
% =========================================================================
function [xrough, yrough] = roughestnew(projection, detector_pixels_xy, thresrough, orderrough, libx, liby, pixel_size, lib_start, lib_end)
    % 物理坐标版粗定位
    nx_pixels = detector_pixels_xy(1);
    ny_pixels = detector_pixels_xy(2);
    x0 = projection(1:nx_pixels);
    y0 = projection(nx_pixels + 1 : nx_pixels + ny_pixels);

    [xdesc0, xindex0] = sort(x0, 'descend');
    [ydesc0, yindex0] = sort(y0, 'descend');

    thres_actual = min(thresrough, nx_pixels);
    weights_x = xdesc0(1:thres_actual).^orderrough;
    weights_y = ydesc0(1:thres_actual).^orderrough;

    xpeak0 = sum(weights_x .* xindex0(1:thres_actual)) / (sum(weights_x) + eps);
    ypeak0 = sum(weights_y .* yindex0(1:thres_actual)) / (sum(weights_y) + eps);
    
    % 转换为物理坐标 (假设中心对齐)
    % 公式: Pos = (Index - (N+1)/2) * Size
    x_phys = (xpeak0 - (nx_pixels + 1) / 2) * pixel_size;
    y_phys = (ypeak0 - (ny_pixels + 1) / 2) * pixel_size;
    
    % 映射到库索引
    lib_step_x = (lib_end - lib_start) / (libx - 1);
    lib_step_y = (lib_end - lib_start) / (liby - 1);
    
    xrough = (x_phys - lib_start) / lib_step_x + 1;
    yrough = (y_phys - lib_start) / lib_step_y + 1;
    
    xrough = max(1, min(libx, round(xrough)));
    yrough = max(1, min(liby, round(yrough)));
end

function [x_interp, y_interp] = interpolate_peak_parabolic_2d(L_matrix, x_idx, y_idx)
    [sub_libx, sub_liby] = size(L_matrix);
    % X
    xm1 = -inf; xp1 = -inf; xc = L_matrix(x_idx, y_idx);
    if x_idx > 1, xm1 = L_matrix(x_idx-1, y_idx); end
    if x_idx < sub_libx, xp1 = L_matrix(x_idx+1, y_idx); end
    denom_x = 2 * (xm1 - 2*xc + xp1);
    if denom_x < 0 && isfinite(xm1) && isfinite(xp1)
        x_interp = x_idx + (xm1 - xp1) / denom_x;
    else
        x_interp = x_idx;
    end
    % Y
    ym1 = -inf; yp1 = -inf; yc = L_matrix(x_idx, y_idx);
    if y_idx > 1, ym1 = L_matrix(x_idx, y_idx-1); end
    if y_idx < sub_liby, yp1 = L_matrix(x_idx, y_idx+1); end
    denom_y = 2 * (ym1 - 2*yc + yp1);
    if denom_y < 0 && isfinite(ym1) && isfinite(yp1)
        y_interp = y_idx + (ym1 - yp1) / denom_y;
    else
        y_interp = y_idx;
    end
end

function [neighbourefer, rangex, rangey] = neighbourdeter_2d(xrough, yrough, referprojection, radius)
    [libx, liby, ~] = size(referprojection);
    
    startx = max(1, xrough - radius);
    endx   = min(libx, xrough + radius);
    rangex = startx : endx;
    
    starty = max(1, yrough - radius);
    endy   = min(liby, yrough + radius);
    rangey = starty : endy;
    
    % 处理 2D 或 3D 库
    if ndims(referprojection) == 2
        neighbourefer = referprojection(rangex, rangey);
    else
        % 如果库是3D的，只取第1层
        neighbourefer = referprojection(rangex, rangey, 1);
    end
end

function L = calc_logL_gaussian(norm_sample_proj, norm_ref_proj)
    pseudolightyield = 1600; 
    norm_ref_proj_safe = norm_ref_proj + eps;
    sample_photons = norm_sample_proj * pseudolightyield;
    ref_photons_safe = norm_ref_proj_safe * pseudolightyield;
    term1 = (sample_photons - ref_photons_safe).^2 ./ (2 * ref_photons_safe);
    term2 = 0.5 * log(2 * pi * ref_photons_safe);
    L = sum(-term1 - term2);
end