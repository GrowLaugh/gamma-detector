% =========================================================================
% === 脚本 2 (批量版): 使用 Anger Logic 批量分析 + 生成汇总图 ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入/输出文件夹 ---

input_data_folder   = 'E:\110work\raw data\calib_data合并';       % 输入文件夹
output_image_folder = 'Anger_Results_ppt';    % 输出图片文件夹
if ~exist(output_image_folder, 'dir'), mkdir(output_image_folder); end

% % --- 文件名格式 ---
% file_basename       = 'filtered_events%04d.mat'; % e.g., calib_events_p0505.mat
file_basename       = 'calib_events_p%04d.mat'; % e.g., calib_events_p0505.mat


% --- !!! 在此处定义所有需要处理的文件编号 !!! ---
file_indices_to_process = [ ...
    503:508, ...  % 第一批
    603:608, ...  % 第二批
    703:708, ...  % 第三批
    803:808  ...  % 第四批
    %     505:506, ...  % 第二批
    % 503,707,909 ...  % 第三批
];

% --- 探测器参数 ---
total_size_mm       = 50;   % 探测器物理尺寸 (50mm x 50mm)
display_resolution  = 256;  % 输出图像的分辨率 (bin数量)

%% 2. 初始化汇总画布
% -------------------------------------------------------------------------
% 预定义直方图边界 (覆盖整个探测器)
axis_limit = total_size_mm / 2; 
bin_edges = linspace(-axis_limit, axis_limit, display_resolution + 1);
% 获取每个bin的中心用于 imagesc
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% *** 初始化主画布 (用于累加所有结果) ***
composite_image = zeros(display_resolution, display_resolution);

%% 3. 主循环：逐个处理文件
% -------------------------------------------------------------------------
fprintf('准备批量处理 %d 个文件...\n', length(file_indices_to_process));
fprintf('============================================================\n');

for file_idx = file_indices_to_process
    
    % --- a. 构建文件名并加载 ---
    current_filename = sprintf(file_basename, file_idx);
    inputFile = fullfile(input_data_folder, current_filename);
    
    if ~exist(inputFile, 'file')
        fprintf('警告: 找不到文件 %s，已跳过。\n', current_filename);
        continue;
    end
    
    fprintf('正在处理: %s ... ', current_filename);
    load(inputFile); % 加载 planeset, pixels_xy
    
    total_events = length(planeset);
    if total_events == 0
        fprintf('无有效事件，跳过。\n');
        continue;
    end
    
    % --- b. 准备 Anger Logic 坐标系 ---
    pixel_size_mm = total_size_mm / pixels_xy(1); 
    pixel_centers_1D = (- (pixels_xy(1)/2 - 0.5) : 1 : (pixels_xy(1)/2 - 0.5)) * pixel_size_mm;
    [x_pos_grid, y_pos_grid] = meshgrid(pixel_centers_1D);
    
    % --- c. 执行 Anger Logic 定位 ---
    all_calc_positions = zeros(total_events, 2);
    for i = 1:total_events
        [cx, cy] = anger_logic(planeset{i}, x_pos_grid, y_pos_grid);
        all_calc_positions(i, :) = [cx, cy];
    end
    
    % --- d. 生成 2D 直方图数据 ---
    % 计算直方图 (Y, X) -> 注意 histcounts2 的输入顺序通常对应矩阵的 (行, 列)
    psf_counts = histcounts2(all_calc_positions(:,2), all_calc_positions(:,1), bin_edges, bin_edges);
    
    % --- e. 累加到汇总图 (归一化) ---
    if max(psf_counts(:)) > 0
        % 归一化当前图，使最高点为1，这样不同强度的源在汇总图中依然清晰可见
        psf_norm = psf_counts / max(psf_counts(:));
        composite_image = composite_image + psf_norm;
    end
    
    % --- f. 绘制并保存单张结果 (可选) ---
    h_fig = figure('Name', 'Anger 2D Result', 'Visible', 'off'); 
    imagesc(bin_centers, bin_centers, psf_counts);
    set(gca, 'YDir', 'normal'); 
    colormap('jet'); colorbar; axis equal tight;
    title(sprintf('Anger Logic: p%04d (N=%d)', file_idx, total_events));
    xlabel('X (mm)'); ylabel('Y (mm)');
    
    % 单张图也加上网格方便查看
    grid on; set(gca, 'GridColor', 'w', 'GridAlpha', 0.4);
    
    output_filename = sprintf('Anger_2D_p%04d.png', file_idx);
    exportgraphics(h_fig, fullfile(output_image_folder, output_filename), 'Resolution', 150);
    close(h_fig);
    
    fprintf('完成\n');
end

fprintf('============================================================\n');

%% 4. 绘制并保存最终汇总图 (含白色网格)
% -------------------------------------------------------------------------
fprintf('正在生成最终汇总图...\n');

h_summary = figure('Name', 'Anger Logic 汇总结果', 'NumberTitle', 'off', 'Position', [100, 100, 800, 700]);

% 绘制累加后的图像
imagesc(bin_centers, bin_centers, composite_image);
set(gca, 'YDir', 'normal'); % 确保Y轴向上为正
colormap('jet'); 
colorbar;
axis equal tight;

% --- 添加白色网格 (5mm 间距) ---
hold on;
grid_ticks = -axis_limit:5:axis_limit; % 从边界开始，每5mm一个刻度
set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks);
grid on;
% 设置网格样式：白色，半透明，实线
set(gca, 'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
set(gca, 'Layer', 'top'); % 确保网格显示在图像上方，不被遮挡
hold off;

% 标题与标签
title(sprintf('Anger Logic 定位结果汇总 (%d 个点源)', length(file_indices_to_process)));
xlabel('X Position (mm)');
ylabel('Y Position (mm)');

% 保存汇总图
summary_filename = 'Anger_2D_Summary_Combined.png';
summary_path = fullfile(output_image_folder, summary_filename);
exportgraphics(h_summary, summary_path, 'Resolution', 300);

fprintf('汇总图已保存至: %s\n', summary_path);
fprintf('所有处理完成！\n');

%% 本地函数
% -------------------------------------------------------------------------
function [x_mm, y_mm] = anger_logic(lightMap, x_grid, y_grid)
%ANGER_LOGIC 加权平均定位
    total_signal = sum(lightMap(:));
    if total_signal == 0
        x_mm = 0; y_mm = 0;
        return;
    end
    sum_x = sum(sum(double(lightMap) .* x_grid));
    sum_y = sum(sum(double(lightMap) .* y_grid));
    
    x_mm = sum_x / total_signal;
    y_mm = sum_y / total_signal;
end