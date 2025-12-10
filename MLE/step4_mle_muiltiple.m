% =========================================================================
% === 脚本 3 (批量版): 使用【混合ML + 3D插值】批量处理、保存多视角图及汇总 ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% % --- 输入文件设置 ---
input_data_folder   ='E:\110work\ExpDataProc\calib_data1';
file_basename       = 'calib_events_p%04d.mat'; % 文件名格式 

% --- !!! 在此处定义所有需要处理的文件编号 !!! ---
file_indices_to_process = [ ...

    503:508, ...  % 第一批
    603:608, ...  % 第二批
    703:708, ...  % 第三批
    802:809,  ...  % 第四批
    901,902,909 ...  % 第三批
];

% --- 输出设置 ---
output_image_folder = 'mle_Results1mm_calib1-new';   % 结果保存文件夹
if ~exist(output_image_folder, 'dir'), mkdir(output_image_folder); end

% --- 参考库设置 ---
libraryFile           ='E:\110work\library生成\projectionslibrary_symmetric_full2Dleijia25.mat';
lib_phys_start        = -25.0; % 参考库起始坐标 (mm)
lib_phys_end          = 25.0;  % 参考库结束坐标 (mm)

grid_interval_mm      = 5;       % 定义白色网格线间距
% --- 显示与网格参数 ---
detector_size_mm      = 50;      % 探测器物理尺寸 (mm)


% --- 定义粗定位的邻域半径 ---
coarse_radius = 5; % 邻域半径 (5 得到一个 11x11 的 *搜索* 窗口)(1mm:r=13)
interp_padding = 1; % 为插值额外计算的边界层

% --- 最终图像显示参数 ---
display_resolution       = 1000;     % 最终图像的分辨率 (像素)
% --- PSF平滑参数 ---
enable_smoothing = true;     % 设置为 true 来启用高斯平滑
smoothing_sigma  = 1;      % 高斯核的标准差 (以bins为单位)

%% 2. 加载参考库 (只加载一次)
% -------------------------------------------------------------------------
fprintf('正在加载参考库: %s ...\n', libraryFile);
if ~exist(libraryFile, 'file')
    error('错误: 找不到参考库文件 %s', libraryFile);
end
load(libraryFile); % 加载 referprojection, x_coords, y_coords, z_coords
[libx, liby, libz] = size(referprojection);
fprintf('参考库加载完成！尺寸: %d x %d x %d\n', libx, liby, libz);

% --- 初始化汇总画布 ---
% 定义物理坐标的边界，确保所有图对齐
axis_limit = detector_size_mm / 2; 
bin_edges = linspace(-axis_limit, axis_limit, display_resolution + 1);
% 获取中心点用于绘图
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% 主画布 (用于累加)
composite_image = zeros(display_resolution, display_resolution);

num_files = length(file_indices_to_process);
fprintf('准备批量处理 %d 个文件...\n', num_files);

% % 启动并行池
% if isempty(gcp('nocreate')), parpool; end

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
    
    % --- b. 执行 ML + 3D 插值定位算法 (并行) ---
    file_calc_indices = zeros(total_events, 3);
    
    parfor i = 1:total_events
        lightMap = planeset{i};
        x_proj = sum(lightMap, 1);
        y_proj = sum(lightMap, 2)';
        combined_proj = [x_proj, y_proj];
        norm_proj = combined_proj / (sum(combined_proj) + eps);
        
        % 1. 粗定位
        try
            % [xrough, yrough] = roughest(norm_proj, pixels_xy, 3, 2, libx, liby);
                    [xrough, yrough] = roughestnew(norm_proj, pixels_xy, 8, 2, ...
                                        libx, liby, 6.0625, ...
                                        lib_phys_start, lib_phys_end);
        catch
            xrough = round(libx/2); yrough = round(liby/2);
        end
        
        % 2. 确定邻域
        [~, rx_s, ry_s] = neighbourdeter(xrough, yrough, referprojection, coarse_radius);
        
        rx_calc_start = max(1, rx_s(1) - interp_padding);
        rx_calc_end   = min(libx, rx_s(end) + interp_padding);
        rx_calc = rx_calc_start:rx_calc_end;
        
        ry_calc_start = max(1, ry_s(1) - interp_padding);
        ry_calc_end   = min(liby, ry_s(end) + interp_padding);
        ry_calc = ry_calc_start:ry_calc_end;
        
        % 3. 3D ML 搜索
        sz_x = length(rx_calc); sz_y = length(ry_calc);
        L_sub = ones(sz_x, sz_y, libz) * -inf;
        
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
        %  % --- 新增步骤 3.5: 对似然矩阵进行空间平滑 ---
        % if enable_ml_smoothing
        %     % 1. 处理 -inf 值：将其替换为矩阵中的有限最小值，避免平滑时产生 NaN
        %     valid_vals = L_sub(isfinite(L_sub));
        %     if ~isempty(valid_vals)
        %         min_val = min(valid_vals);
        %         L_sub(L_sub == -inf) = min_val - 10; % 设置一个比最小值还低的底板
        % 
        %         % 2. 执行 3D 高斯平滑
        %         % 如果没有 Image Processing Toolbox，可以使用 fallback 方法
        %         try
        %             L_sub = imgaussfilt3(L_sub, ml_smoothing_sigma);
        %         catch
        %             % 降级处理：只做 2D 平滑 (逐层 Z)
        %             for z_layer = 1:libz
        %                 L_sub(:,:,z_layer) = imgaussfilt(L_sub(:,:,z_layer), ml_smoothing_sigma);
        %             end
        %         end
        %     end
        % end

        % 4. 寻找峰值
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
            % 5. 抛物线插值
            [ix, iy, iz] = interpolate_peak_parabolic(L_sub, sx, sy, sz);
            file_calc_indices(i, :) = [ix + rx_calc_start - 1, ...
                                       iy + ry_calc_start - 1, ...
                                       iz];
        else
            file_calc_indices(i, :) = [NaN, NaN, NaN];
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
    dz = z_coords(2)-z_coords(1); z0 = z_coords(1);
    
    xpo_mm = (file_calc_indices(:,1) - 1) * dx + x0;
    ypo_mm = (file_calc_indices(:,2) - 1) * dy + y0;
    zpo_mm = (file_calc_indices(:,3) - 1) * dz + z0;
    
    % --- d. 生成直方图数据 ---
    h_obj = histogram2(xpo_mm, ypo_mm, bin_edges, bin_edges, 'Visible', 'off');
    current_img = h_obj.Values'; % 转置以匹配 imagesc
    

    
    % ================= 保存图片部分 =================
    
    % --- 1. 保存 2D 结果图 (带网格) ---
    fig_2d = figure('Name', '2D Result', 'Visible', 'off');
    imagesc(bin_centers, bin_centers, current_img);
    set(gca, 'YDir', 'normal');
    colormap('jet'); colorbar; axis equal tight;
    title(sprintf('2D MLE: p%04d (N=%d)', file_idx, length(xpo_mm)));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % 添加白色网格 (5mm)
    hold on;
    grid_ticks = -25:5:25;
    set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks);
    grid on;
    set(gca, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
    set(gca, 'Layer', 'top');
    hold off;
    
    exportgraphics(fig_2d, fullfile(output_image_folder, sprintf('MLE_2D_p%04d.png', file_idx)));
    close(fig_2d);
    
    % --- 2. 保存 3D 结果图 ---
    fig_3d = figure('Name', '3D Result', 'Visible', 'off');
    histogram2(xpo_mm, ypo_mm, bin_edges, bin_edges, 'FaceColor', 'flat', 'ShowEmptyBins', 'off');
    view(45, 30); 
    colormap('jet'); colorbar;
    title(sprintf('3D MLE: p%04d', file_idx));
    xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Counts');
    
    exportgraphics(fig_3d, fullfile(output_image_folder, sprintf('MLE_3D_p%04d.png', file_idx)));
    close(fig_3d);
    
    % --- 3. 保存 Z (DOI) 分布图 ---
    fig_z = figure('Name', 'Z Distribution', 'Visible', 'off');
    histogram(zpo_mm, libz*4); % 增加些bins看细节
    grid on;
    title(sprintf('Depth of Interaction (Z): p%04d', file_idx));
    xlabel('Z Position (mm)'); ylabel('Counts');
    
    exportgraphics(fig_z, fullfile(output_image_folder, sprintf('MLE_Z_p%04d.png', file_idx)));
    close(fig_z);
    
    % --- e. 累加到汇总图 ---
    if max(current_img(:)) > 0
        % 归一化后累加，使每个点源亮度一致
        composite_image = composite_image + (current_img / max(current_img(:)));
    end
    
    fprintf('处理完毕，图片已保存。\n');
end

fprintf('============================================================\n');

%% 4. 生成并保存最终汇总图
% -------------------------------------------------------------------------
fprintf('正在生成汇总图...\n');
h_sum = figure('Name', 'Summary', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);

imagesc(bin_centers, bin_centers, composite_image);
set(gca, 'YDir', 'normal');
colormap('jet'); colorbar; axis equal tight;

% --- 添加白色网格 (5mm) ---
hold on;
grid_ticks = -25:5:25; 
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
function [x_interp, y_interp, z_interp] = interpolate_peak_parabolic(L_matrix, x_idx, y_idx, z_idx)
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
    pseudolightyield = 8000; 
    norm_ref_proj_safe = norm_ref_proj + eps;
    sample_photons = norm_sample_proj * pseudolightyield;
    ref_photons_safe = norm_ref_proj_safe * pseudolightyield;
    term1 = (sample_photons - ref_photons_safe).^2 ./ (2 * ref_photons_safe);
    term2 = 0.5 * log(2 * pi * ref_photons_safe);
    L = sum(-term1 - term2);
end