% =========================================================================
% === 脚本 1 (批量版): 批量处理多个.bin文件 (全程手动筛选 & 独立保存) ===
% =========================================================================
clear;
clc;
close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入文件设置 ---
input_data_folder   = 'E:\110work\raw data\raw data1';       % 存放 .bin 文件的文件夹

% 文件名格式, e.g., p0505.bin, p0708.bin etc.
file_basename       = 'p%04d.bin'; 

% --- !!! 关键：在此处定义所有需要处理的文件编号 !!! ---
file_indices_to_process = [ ...
    503:508, ...  % 第一批
    603:608, ...  % 第二批
    703:708, ...  % 第三批
    801:809,  ...  % 第四批
    901,902,909 ...  % 第三批
];

% --- 输出文件设置 ---
output_folder       = 'calib_data1';       % 保存 .mat 文件的文件夹
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

% --- 能谱图片保存设置 ---
energy_spec_folder  = 'energy_spec1'; % 保存能谱图的文件夹
if ~exist(energy_spec_folder, 'dir'), mkdir(energy_spec_folder); end

% --- 数据处理参数 ---
nChannels       = 64;
pixels_xy       = [8, 8];
adc_conversion_factor = 0.0012207;

%% 2. 加载ADC校准数据
% -------------------------------------------------------------------------
use_calibration = false;
if exist('Offsetdata.mat', 'file') && exist('Gaindata.mat', 'file')
    load('Offsetdata.mat', 'offset_value'); 
    load('Gaindata.mat', 'gain_value');
    use_calibration = true;
    fprintf('成功加载 ADC 校准文件 (Offset & Gain)。\n');
else
    warning('未找到 Offsetdata.mat 或 Gaindata.mat，将跳过 ADC 校准步骤！');
end

fprintf('准备开始批量处理 %d 个文件...\n', length(file_indices_to_process));
fprintf('注意：本模式需要您对【每一个】文件手动点击选择能量窗。\n');
fprintf('============================================================\n');

%% 3. 主循环：逐个处理每个.bin文件
% -------------------------------------------------------------------------
for file_idx = file_indices_to_process
    
    % --- a. 构建并加载当前文件 ---
    current_bin_filename = sprintf(file_basename, file_idx);
    inputFile = fullfile(input_data_folder, current_bin_filename);
    
    if ~exist(inputFile, 'file')
        fprintf('警告: 找不到文件 %s，已跳过。\n', inputFile);
        continue;
    end
    
    fprintf('------------------------------------------------------------\n');
    fprintf('正在处理文件: %s ...\n', current_bin_filename);
    
    % --- b. 读取数据 ---
    fid = fopen(inputFile);
    if fid < 0, continue; end
    rawData = fread(fid, [nChannels, Inf], '64*int16', 8);
    fclose(fid);
    
    rawData = rawData * adc_conversion_factor;
    n_raw_events = size(rawData, 2);

    % --- c. 应用 ADC 校准 ---
    if use_calibration
        data_calib = zeros(size(rawData));
        for k = 1:n_raw_events
            data_calib(:, k) = (rawData(:, k) - offset_value) ./ gain_value;
        end
    else
        data_calib = rawData;
    end
    
    energy = sum(data_calib, 1);
    
    % --- d. 能谱绘制与手动筛选 (核心修改) ---
    
    % 1. 创建可见的图形窗口
    h_fig = figure('Name', sprintf('Energy Spectrum - %s', current_bin_filename), 'Visible', 'on');
    
    % 2. 绘制直方图
    c1 = 0:0.2:100; 
    c = hist(energy, c1);
    bar(c1, c, 'histc');
    xlabel('能量 (单位)'); ylabel('计数');
    title(sprintf('Energy Spectrum: %s (Total: %d) - 请点击两次设定窗口', current_bin_filename, n_raw_events));
    grid on; hold on;
    
    % 3. 强制刷新并等待交互
    figure(h_fig); % 将焦点给到这个窗口
    drawnow;       % 强制立即渲染
    pause(0.2);    % 防止白屏，确保渲染完成
    
    fprintf('  >>> 请在图中点击两次：1.左边界(下限)  2.右边界(上限)...\n');
    
    % 4. 获取鼠标点击
    [x_clicks, ~] = ginput(2);
    eLw = min(x_clicks);
    eHw = max(x_clicks);
    
    fprintf('  能量窗已设定为: [%.2f, %.2f]\n', eLw, eHw);
    
    % 5. 在图上绘制阈值线
    xline(eLw, 'g--', 'LineWidth', 2);
    xline(eHw, 'g--', 'LineWidth', 2);
    
    % 6. 保存能谱图
    spec_img_name = strrep(current_bin_filename, '.bin', '.png');
    spec_img_path = fullfile(energy_spec_folder, ['spec_' spec_img_name]);
    exportgraphics(h_fig, spec_img_path);
    
    % 7. 关闭窗口
    close(h_fig);
    
    % --- e. 应用筛选 ---
    valid_idx = find(energy >= eLw & energy <= eHw);
    data_valid = data_calib(:, valid_idx);
    
    num_events_in_file = size(data_valid, 2);
    fprintf('  > 原始: %d -> 有效: %d (%.1f%%) | 能谱已存至 %s\n', ...
        n_raw_events, num_events_in_file, 100*num_events_in_file/n_raw_events, energy_spec_folder);

    % --- f. 提取事件并保存为独立文件 ---
    if num_events_in_file > 0
        % 1. 转换为 planeset
        planeset = cell(1, num_events_in_file);
        for j = 1:num_events_in_file
            lightMap = reshape(data_valid(:, j), pixels_xy(1), pixels_xy(2));
            planeset{j} = lightMap;
        end
        
        % 2. 构建输出文件名 (calib_events_pXXXX.mat)
        current_output_filename = sprintf('calib_events_p%04d.mat', file_idx);
        current_output_path = fullfile(output_folder, current_output_filename);
        
        % 3. 保存单个文件 (包含该文件特定的阈值 eLw, eHw)
        save(current_output_path, 'planeset', 'pixels_xy', 'eLw', 'eHw');
        fprintf('  > 已保存有效事件到: %s\n', current_output_path);
    else
        fprintf('  > 警告: 此文件没有有效事件，跳过保存。\n');
    end
end

fprintf('============================================================\n');
fprintf('所有文件处理完毕！\n');
fprintf('  - 结果保存在文件夹: %s\n', output_folder);
fprintf('  - 能谱图保存在文件夹: %s\n', energy_spec_folder);