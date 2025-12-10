clear all;
close all;

%% 1. 读取数据
data_folder='E:\110work\raw data\raw data1';
input_filename='p0505.bin';
inputFile             = fullfile(data_folder,input_filename);

fid = fopen(inputFile,'rb');

% 把源文件的N*64变为64*N,这样有64行，每一行代表某位置处的光子数记录    [64个数据] [8字节头] [64个数据] ...。
% data1 = fread(fid,[64 ,Inf],"int16=>double",8); % 每次读64通道(64个数值)，再跳过8字节，再读下一个64通道，64*int16 表示每次读取64个16位整数
data1 = fread(fid,[64 ,Inf],'64*int16=>double',8); 

fclose(fid);

% 将读取的数据转换成能量单位
data1 = data1 * 0.0012207;

% 沿着行求和，每一行为一个探测器通道的信号强度
energy = sum(data1,1);  % 对每一列求和，得到一个能量向量 energy，其长度为事件的数量，每个元素是单事件的能量。

%% 绘制能量直方图，选择右侧峰的起始作为能量阈值
c1=0:0.1:20;
c=hist(energy,c1);
% 绘制直方图
figure;
bar(c1, c, 'histc');  % 使用 'histc' 来生成一个带有条形的直方图
xlabel('能量 (单位)');
ylabel('频数');
title('能量直方图');
grid on;  % 打开网格

%% 2. 能量筛选
idx_energy = find(energy > 11);   % 直接找到能量大于11的事件
data2 = data1(:, idx_energy);%data2 仍然是 64 × N 矩阵，每一列对应一个事件


%% 4. 整理为 8×8×N 事件矩阵
nKeep = size(data2,2);
data5 = reshape(data2,[8,8,nKeep]);

%% 5. 可视化（空格翻页查看事件）
figure('Name', '矩阵可视化', 'NumberTitle', 'off', ...
       'KeyPressFcn', @(src,event) setappdata(gcf, 'key', event.Key));

currentIndex = 1; % 显示第一组数据

while ishandle(gcf) % 当窗口存在时循环
    key = getappdata(gcf, 'key');
    if strcmp(key, 'space')
        currentIndex = mod(currentIndex, nKeep) + 1;
        data = squeeze(data5(:,:,currentIndex));

        % 清除当前图形并绘制新图
        clf;
        bar3(data);

        % 添加基本标签
        title(sprintf('矩阵组号: %d/%d', currentIndex, nKeep), 'FontSize', 12);
        xlabel('列索引');
        ylabel('行索引');
        zlabel('数值');

        % 设置一致的Z轴范围
        zlim([min(data5(:)), 5]);

        % 添加操作提示
        annotation('textbox', [0.05, 0.01, 0.9, 0.04], ...
                   'String', '按空格键显示下一组', ...
                   'EdgeColor', 'none', ...
                   'FontSize', 10, ...
                   'HorizontalAlignment', 'center');

        drawnow; % 更新显示
        setappdata(gcf, 'key', ''); % 重置按键状态
    end
    pause(0.05); % 短暂停避免CPU过载
end





% %% 3. 位置筛选
% nEvents = size(data2,2);%获取data2的列数，即有效事件的总数
% keepIdx = false(1,nEvents);   % 逻辑索引，标记符合条件的事件
% 
% for i = 1:nEvents %遍历所有事件
%     data3 = reshape(data2(:,i),8,8); %将每个事件的数据重新塑造成 8x8 的矩阵
%     total = sum(data3(:)); %计算该事件的总能量
% 
%     % 计算重心
%     row_weights = sum(data3, 2);%对每一行求和
%     peak_row = (1:8) * row_weights / total;%行号 × 对应的信号强度，再把所有行加起来，再除以总能量
% 
%     col_weights = sum(data3, 1);
%     peak_col = (1:8) * col_weights' / total;
% end
% 
% 

% =========================================================================
% === 脚本：预处理伽马相机实验数据 (.bin) 版本 ===
% =========================================================================
% clear; clc; close all;
% 
% %% 1. 参数设置
% % -------------------------------------------------------------------------
% inputFile       = 'p0706.bin';   % 输入BIN文件
% nChannels       = 64;                 % 探测器通道数 (8x8)
% photonThreshold = 15;                 % 能量阈值 (需根据实验调节)
roi_center      = [4.5, 4.5];         % 感兴趣区域中心 (row, col)
roi_half_range  = 5.0;                % ROI半径 (默认=1 → 中间 4x4)
% outputScintiFile   = 'scintisetsamples.mat';
% outputPlanesetFile = 'planesetsamples.mat';
% 
% %% 2. 读取BIN文件
% % -------------------------------------------------------------------------
% fprintf('正在读取实验数据: %s ...\n', inputFile);
% fid = fopen(inputFile, 'rb');
% if fid < 0
%     error('无法打开文件: %s', inputFile);
% end
% rawData = fread(fid, [nChannels, Inf], 'int16=>double', 8); % 每次64通道
% fclose(fid);
% 
% % ADC转换系数 (根据实验设置)
% rawData = rawData * 0.0012207;  
% 
% fprintf('文件读取完成，共 %d 个采样事件。\n', size(rawData, 2));
% 
% %% 3. 能量计算和阈值筛选
% % -------------------------------------------------------------------------
% energy = sum(rawData, 1);   % 每次事件的能量
% valid_idx = find(energy > photonThreshold);
% data_valid = rawData(:, valid_idx);
% fprintf('能量筛选后剩余 %d 个事件。\n', size(data_valid, 2));
% 
%% 4. 重心计算 + ROI筛选
% -------------------------------------------------------------------------
nEvents = size(data_valid, 2);
scinti  = {};   % 存储“源点位置”（用重心代替）
planeset = {};  % 存储每次事件的 8x8 光分布矩阵

for i = 1:nEvents
    data3 = reshape(data_valid(:, i), 8, 8); % 8x8 矩阵
    total = sum(data3(:));
    if total == 0, continue; end

    % 行/列加权平均
    row_weights = sum(data3, 2);  
    peak_row = (1:8) * row_weights / total;
    col_weights = sum(data3, 1);  
    peak_col = (1:8) * col_weights' / total;

    % ROI筛选
    if (peak_row > roi_center(1)-roi_half_range) && (peak_row < roi_center(1)+roi_half_range) && ...
       (peak_col > roi_center(2)-roi_half_range) && (peak_col < roi_center(2)+roi_half_range)

        scinti{end+1}  = [peak_row, peak_col]; % 源点 (重心坐标)
        planeset{end+1} = data3;               % 光分布 (8x8)
    end
end

fprintf('ROI筛选后剩余 %d 个有效事件。\n', length(scinti));

%% 5. 保存结果
% -------------------------------------------------------------------------
save(outputScintiFile, 'scinti');
save(outputPlanesetFile, 'planeset');
fprintf('处理完成！结果已保存到:\n - %s\n - %s\n', ...
    outputScintiFile, outputPlanesetFile);

