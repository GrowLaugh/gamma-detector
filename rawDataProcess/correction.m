% =========================================================================
% === 脚本 0: 将原始 Gain/Offset 二进制文件转化为校准用 .mat 文件 ===
% =========================================================================
clear; clc; close all;

%% 1. 参数设置
% -------------------------------------------------------------------------
% --- 输入路径设置 (请修改为您存放 Gain.bin 和 Offset.bin 的文件夹) ---
input_data_folder = 'E:\110work\raw data\raw data2'; 

% --- 文件名 ---
offset_bin_file = 'offset.bin';
gain_bin_file   = 'GAIN1.bin';

% --- 输出文件名 ---
output_offset_mat = 'Offsetdata.mat';
output_gain_mat   = 'Gaindata.mat';

% --- 数据参数 ---
nChannels = 64;
% ADC转换系数 (必须与处理实验数据时保持一致!)
adc_conversion_factor = 0.0012207; 

%% 2. 处理 Offset 文件
% -------------------------------------------------------------------------
offset_path = fullfile(input_data_folder, offset_bin_file);
fprintf('------------------------------------------------------------\n');
fprintf('正在处理 Offset 文件: %s ...\n', offset_path);

if exist(offset_path, 'file')
    % 读取二进制数据
    fid = fopen(offset_path, 'rb');
    % 格式: 64通道 int16, 每次跳过8字节头
    raw_offset_data = fread(fid, [nChannels, Inf], '64*int16', 8);
    fclose(fid);
    
    n_offset_events = size(raw_offset_data, 2);
    fprintf('  读取到 %d 个 Offset 采样事件。\n', n_offset_events);
    
    % 转换为物理单位 (与主程序保持一致)
    raw_offset_data = raw_offset_data * adc_conversion_factor;
    
    % 计算每个通道的平均值作为 Offset
    % offset_value 应该是一个 64x1 的向量
    offset_value = mean(raw_offset_data, 2);
    
    % 保存 Offsetdata.mat
    save(output_offset_mat, 'offset_value');
    fprintf('  >>> Offset 校准文件已保存: %s\n', output_offset_mat);
    
    % 可视化 Offset
    figure('Name', 'Offset Calibration');
    bar(offset_value);
    title('Offset per Channel');
    xlabel('Channel Index'); ylabel('Offset Value');
    grid on;
else
    error('错误: 找不到 Offset 文件: %s', offset_path);
end

%% 3. 处理 Gain 文件
% -------------------------------------------------------------------------
gain_path = fullfile(input_data_folder, gain_bin_file);
fprintf('------------------------------------------------------------\n');
fprintf('正在处理 Gain 文件: %s ...\n', gain_path);

if exist(gain_path, 'file')
    % 读取二进制数据
    fid = fopen(gain_path);
    raw_gain_data = fread(fid, [nChannels, Inf], '64*int16', 8);
    fclose(fid);
    
    n_gain_events = size(raw_gain_data, 2);
    fprintf('  读取到 %d 个 Gain 采样事件。\n', n_gain_events);
    
    % 转换为物理单位
    raw_gain_data = raw_gain_data * adc_conversion_factor;
    
    % 计算每个通道的原始平均值
    gain_raw_mean = mean(raw_gain_data, 2);
    
    % --- 关键步骤: 计算净增益 (Net Gain) ---
    % 1. 减去 Offset (得到纯净的信号幅度)
    gain_net = gain_raw_mean - offset_value;
    
    % 2. 处理死像素或异常值 (防止除以0)
    gain_net(gain_net <= 0) = mean(gain_net(gain_net > 0)); 
    
    % 3. 归一化 (Normalization)
    % 计算所有通道的平均响应强度
    global_mean_response = mean(gain_net);
    
    % 计算相对增益系数 (Relative Gain Factor)
    % 这样处理后，gain_value 的均值为 1。
    % 灵敏度高的通道 gain > 1，灵敏度低的通道 gain < 1。
    gain_value = gain_net / global_mean_response;
    
    % 保存 Gaindata.mat
    % 同时保存 absolute_gain 以备不时之需
    gain_absolute = gain_net; 
    save(output_gain_mat, 'gain_value', 'gain_absolute', 'global_mean_response');
    fprintf('  >>> Gain 校准文件已保存: %s\n', output_gain_mat);
    
    % 可视化 Gain
    figure('Name', 'Gain Calibration');
    subplot(2,1,1);
    bar(gain_net); title('Net Signal Amplitude (Mean - Offset)'); grid on;
    subplot(2,1,2);
    bar(gain_value); title('Normalized Gain Factors (Mean=1)'); grid on;
    yline(1, 'r--', 'LineWidth', 2);
else
    error('错误: 找不到 Gain 文件: %s', gain_path);
end

fprintf('------------------------------------------------------------\n');
fprintf('校准文件生成完毕！请确保将生成的 .mat 文件放在主脚本的工作目录下。\n');