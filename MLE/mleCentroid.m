% =========================================================================
% === 脚本 3 (批量版): 使用【2D MLE + 质心法】批量处理 (抗网格化版) ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入文件设置 ---
input_data_folder   = 'E:\110work\ExpDataProc\calib_data1';
file_basename       = 'calib_events_p%04d.mat'; 

file_indices_to_process = [ ...
   
    503:508, ... 
    603:608, ... 
    703:708, ... 
    801:809, ...  
    901, 902, 909 ...
];

% --- 输出设置 ---
output_image_folder = 'mle_Results1mm_2D_Centroid-pixel256-48.5';   
if ~exist(output_image_folder, 'dir'), mkdir(output_image_folder); end

% --- 参考库设置 ---
libraryFile = 'E:\110work\library生成\projectionslibrary_symmetric_full2D25_Synthetic.mat';

% --- 物理坐标参数 ---
lib_phys_start        = -25.0; 
lib_phys_end          = 25.0;  
detector_size_mm      = 48.5;    
pixel_size_mm         = 6.0625; 

% --- 算法参数 ---
coarse_radius    = 5;  
interp_padding   = 2;  % (重要) 质心法建议边界稍大，例如2或3

% --- !!! 核心改进：使用质心法抗网格化 !!! ---
% 'centroid': 质心法 (推荐，解决5x5方块问题)
% 'parabolic': 抛物线法 (旧方法)
positioning_method    = 'centroid';

% 温度系数 (Temperature, T) - 仅用于 centroid
% 用来"软化"极度尖锐的似然分布。
% 如果结果依然网格化(5x5方块)，增大 T (e.g. 50, 100)。
% 如果结果太糊，减小 T (e.g. 5, 10)。
centroid_temperature  = 6.0; 

% 似然矩阵平滑 (辅助抗噪)
enable_ml_smoothing   = false;    
ml_smoothing_sigma    = 0.8;     

% --- 显示与网格参数 ---
display_resolution    = 256;    
grid_interval_mm      = 5;       

% --- 结果图显示平滑 ---
enable_result_smoothing = false;  
result_smoothing_sigma  = 2.0; % 稍微增大一点显示平滑，让1000分辨率下的点更好看   

%% 2. 加载参考库
% -------------------------------------------------------------------------
fprintf('正在加载参考库: %s ...\n', libraryFile);
if ~exist(libraryFile, 'file'), error('错误: 找不到参考库文件'); end
load(libraryFile); 

[libx, liby, libz] = size(referprojection);
if libz > 1, libz = 1; end % 强制2D
fprintf('参考库加载成功！尺寸: %d x %d (2D)\n', libx, liby);

% --- 初始化汇总画布 ---
axis_limit = detector_size_mm / 2; 
bin_edges = linspace(-axis_limit, axis_limit, display_resolution + 1);
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
composite_image = zeros(display_resolution, display_resolution); 

num_files = length(file_indices_to_process);
fprintf('准备批量处理 %d 个文件...\n', num_files);

if isempty(gcp('nocreate')), parpool; end

%% 3. 主循环
% -------------------------------------------------------------------------
for k = 1:num_files
    file_idx = file_indices_to_process(k);
    current_filename = sprintf(file_basename, file_idx);
    inputFile = fullfile(input_data_folder, current_filename);
    
    if ~exist(inputFile, 'file')
        fprintf('警告: 找不到文件 %s\n', current_filename); continue;
    end
    
    fprintf('------------------------------------------------------------\n');
    fprintf('正在处理 (%d/%d): %s ... ', k, num_files, current_filename);
    
    data_struct = load(inputFile);
    planeset = data_struct.planeset;
    pixels_xy = data_struct.pixels_xy;
    total_events = length(planeset);
    
    if total_events == 0, fprintf('无有效事件。\n'); continue; end
    
    % --- 执行 2D ML 定位 ---
    file_calc_indices = zeros(total_events, 2);
    
    parfor i = 1:total_events
        lightMap = planeset{i};
        x_proj = sum(lightMap, 1);
        y_proj = sum(lightMap, 2)';
        
        combined_proj = [x_proj, y_proj];
        norm_proj = combined_proj / (sum(combined_proj) + eps);
        
        % 1. 粗定位
        try
            [xrough, yrough] = roughestnew(norm_proj, pixels_xy, 8, 2, ...
                                           libx, liby, pixel_size_mm, ...
                                           lib_phys_start, lib_phys_end);
        catch
            xrough = round(libx/2); yrough = round(liby/2);
        end
        
        % 2. 确定邻域
        [~, rx_s, ry_s] = neighbourdeter_2d(xrough, yrough, referprojection, coarse_radius);
        
        % 扩大范围
        rx_calc_start = max(1, rx_s(1) - interp_padding);
        rx_calc_end   = min(libx, rx_s(end) + interp_padding);
        rx_calc = rx_calc_start:rx_calc_end;
        
        ry_calc_start = max(1, ry_s(1) - interp_padding);
        ry_calc_end   = min(liby, ry_s(end) + interp_padding);
        ry_calc = ry_calc_start:ry_calc_end;
        
        % 3. 2D ML 搜索
        sz_x = length(rx_calc); 
        sz_y = length(ry_calc);
        L_sub = ones(sz_x, sz_y) * -inf; 
        
        for xi = 1:sz_x
            for yi = 1:sz_y
                xg = rx_calc(xi); yg = ry_calc(yi);
                if ndims(referprojection) == 2
                    ref = referprojection{xg, yg};
                else
                    ref = referprojection{xg, yg, 1};
                end
                
                if ~isempty(ref)
                    L_sub(xi, yi) = calc_logL_gaussian(norm_proj, ref);
                end
            end
        end
        
        % 4. 平滑处理
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
        
        % 5. 精确定位 (算法切换)
        if strcmp(positioning_method, 'centroid')
            % === 推荐：质心法 (解决5x5方块) ===
            [ix, iy] = calculate_mle_centroid_2d(L_sub, centroid_temperature);
        else
            % === 旧方法：抛物线 ===
            [max_L, lin_idx] = max(L_sub(:));
            if max_L > -inf
                [sx, sy] = ind2sub(size(L_sub), lin_idx);
                [ix, iy] = interpolate_peak_parabolic_2d(L_sub, sx, sy);
            else
                ix=NaN; iy=NaN;
            end
        end
        
        % 6. 转换索引
        if ~isnan(ix)
            file_calc_indices(i, :) = [ix + rx_calc_start - 1, ...
                                       iy + ry_calc_start - 1];
        else
            file_calc_indices(i, :) = [NaN, NaN];
        end
    end
    
    valid_mask = ~isnan(file_calc_indices(:,1));
    file_calc_indices = file_calc_indices(valid_mask, :);
    
    if isempty(file_calc_indices), fprintf('定位失败。\n'); continue; end
    
    % --- c. 转换为物理坐标 ---
    dx = (lib_phys_end - lib_phys_start) / (libx - 1); x0 = lib_phys_start;
    dy = (lib_phys_end - lib_phys_start) / (liby - 1); y0 = lib_phys_start;
    
    xpo_mm = (file_calc_indices(:,1) - 1) * dx + x0;
    ypo_mm = (file_calc_indices(:,2) - 1) * dy + y0;
    
    % --- d. 生成直方图数据 (histcounts2) ---
    raw_counts = histcounts2(ypo_mm, xpo_mm, bin_edges, bin_edges);
    current_img = raw_counts; 
    
    if enable_result_smoothing
        current_img = imgaussfilt(current_img, result_smoothing_sigma);
    end
    
    % --- 保存图片 ---
    
    % 2D
    fig_2d = figure('Name', '2D Result', 'Visible', 'off');
    imagesc(bin_centers, bin_centers, current_img);
    set(gca, 'YDir', 'normal'); colormap('jet'); colorbar; axis equal tight;
    title(sprintf('2D MLE (%s): p%04d', positioning_method, file_idx));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    hold on;
    grid_ticks = -axis_limit:grid_interval_mm:axis_limit;
    set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
    grid on; set(gca, 'Layer', 'top'); hold off;
    
    exportgraphics(fig_2d, fullfile(output_image_folder, sprintf('MLE_2D_p%04d.png', file_idx)));
    close(fig_2d);
    
    % 3D
    fig_3d = figure('Name', '3D Result', 'Visible', 'off');
    histogram2(xpo_mm, ypo_mm, bin_edges, bin_edges, 'FaceColor', 'flat', 'ShowEmptyBins', 'off');
    view(45, 30); colormap('jet'); colorbar;
    title(sprintf('3D MLE: p%04d', file_idx));
    exportgraphics(fig_3d, fullfile(output_image_folder, sprintf('MLE_3D_p%04d.png', file_idx)));
    close(fig_3d);
    
    % 汇总
    if max(current_img(:)) > 0
        composite_image = composite_image + (current_img / max(current_img(:)));
    end
    %     % 归一化汇总图，确保 Colorbar 最大值为 1
    % if max(composite_image(:)) > 0
    %     composite_image = composite_image / max(composite_image(:));
    % end


    fprintf('完成。\n');
end

fprintf('============================================================\n');
fprintf('生成最终汇总图...\n');
h_sum = figure('Name', 'Summary', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);
imagesc(bin_centers, bin_centers, composite_image);
set(gca, 'YDir', 'normal'); colormap('jet'); colorbar; axis equal tight;

hold on;
grid_ticks = -axis_limit:grid_interval_mm:axis_limit;
set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
grid on; set(gca, 'Layer', 'top'); hold off;

title(sprintf('MLE汇总 (%s, T=%.1f)', positioning_method, centroid_temperature));
exportgraphics(h_sum, fullfile(output_image_folder, 'MLE_Summary_Combined.png'), 'Resolution', 300);
fprintf('任务完成！\n');

%% 本地函数
% =========================================================================
function [xrough, yrough] = roughestnew(projection, detector_pixels_xy, thresrough, orderrough, libx, liby, pixel_size, lib_start, lib_end)
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
    
    x_phys = (xpeak0 - (nx_pixels + 1) / 2) * pixel_size;
    y_phys = (ypeak0 - (ny_pixels + 1) / 2) * pixel_size;
    
    lib_step_x = (lib_end - lib_start) / (libx - 1);
    lib_step_y = (lib_end - lib_start) / (liby - 1);
    
    xrough = (x_phys - lib_start) / lib_step_x + 1;
    yrough = (y_phys - lib_start) / lib_step_y + 1;
    
    xrough = max(1, min(libx, round(xrough)));
    yrough = max(1, min(liby, round(yrough)));
end

function [x_cent, y_cent] = calculate_mle_centroid_2d(L_matrix, T)
% 使用带温度系数的质心法
    L_max = max(L_matrix(:));
    if isinf(L_max), x_cent=NaN; y_cent=NaN; return; end
    
    % 核心公式: W = exp( (L - Lmax) / T )
    W = exp((L_matrix - L_max) / T);
    W_total = sum(W(:));
    
    if W_total == 0, x_cent=NaN; y_cent=NaN; return; end
    
    [sz_x, sz_y] = size(W);
    [X_grid, Y_grid] = ndgrid(1:sz_x, 1:sz_y);
    
    x_cent = sum(W(:) .* X_grid(:)) / W_total;
    y_cent = sum(W(:) .* Y_grid(:)) / W_total;
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
    startx = max(1, xrough - radius); endx = min(libx, xrough + radius); rangex = startx : endx;
    starty = max(1, yrough - radius); endy = min(liby, yrough + radius); rangey = starty : endy;
    if ndims(referprojection) == 2
        neighbourefer = referprojection(rangex, rangey);
    else
        neighbourefer = referprojection(rangex, rangey, 1);
    end
end

% function L = calc_logL_gaussian(norm_sample_proj, norm_ref_proj)
%     pseudolightyield = 8000; 
%     norm_ref_proj_safe = norm_ref_proj + eps;
%     sample_photons = norm_sample_proj * pseudolightyield;
%     ref_photons_safe = norm_ref_proj_safe * pseudolightyield;
%     term1 = (sample_photons - ref_photons_safe).^2 ./ (2 * ref_photons_safe);
%     term2 = 0.5 * log(2 * pi * ref_photons_safe);
%     L = sum(-term1 - term2);
% end
function L = calc_logL_gaussian(norm_sample_proj, norm_ref_proj)
% =========================================================================
%  *** ML 核心函数 (高斯近似) - NaI(Tl) @ 140keV 版本 ***
% =========================================================================
    
    % 估算依据: 
    % 1. NaI(Tl) 光产额 ~38,000/MeV
    % 2. 能量 140keV (0.14 MeV) -> 产生约 5320 光子
    % 3. 探测效率 ~30% -> 最终约 1600 光子
    pseudolightyield = 1600; 
    
    % (以下逻辑不变)
    norm_ref_proj_safe = norm_ref_proj + eps;
    sample_photons = norm_sample_proj * pseudolightyield;
    ref_photons_safe = norm_ref_proj_safe * pseudolightyield;
    
    term1 = (sample_photons - ref_photons_safe).^2 ./ (2 * ref_photons_safe);
    term2 = 0.5 * log(2 * pi * ref_photons_safe);
    
    L = sum(-term1 - term2);
end