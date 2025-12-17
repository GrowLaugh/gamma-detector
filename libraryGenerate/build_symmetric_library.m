% =========================================================================
% === 脚本：读取1/4参考库文件，并将其对称扩展为全尺寸库 (修正轴线覆盖问题) ===
% =========================================================================
clear;
clc;
close all;

%% 1. 定义输入和输出文件
% -------------------------------------------------------------------------
inputFile = 'MyReferenceLibrary3mm.mat'; 
outputFile = 'new48.5mmstep3mm.mat';

%% 2. 加载已有的1/4参考库
% -------------------------------------------------------------------------
fprintf('正在加载已有的1/4参考库: %s ...\n', inputFile);
if ~exist(inputFile, 'file')
    error('错误：找不到输入文件！请确保 %s 文件与此脚本在同一个目录下。', inputFile);
end
load(inputFile);
fprintf('文件加载成功！\n');

%% 3. 初始化全尺寸库
% -------------------------------------------------------------------------
referprojection_q1 = referprojection;
lightMapLibrary_q1 = lightMapLibrary;
x_coords_quad1 = x_coords;
y_coords_quad1 = y_coords;

% 获取第一象限库的维度
[nx_quad, ny_quad, nz] = size(referprojection_q1);

% 计算完整库的维度
nx_full = nx_quad * 2 - 1;
ny_full = ny_quad * 2 - 1;
% nx_full = nx_quad * 2 ;
% ny_full = ny_quad * 2 ;%针对起始坐标不是0的情况



% 创建代表完整物理范围的新坐标轴
x_coords_full = [-fliplr(x_coords_quad1(2:end)), x_coords_quad1];
y_coords_full = [-fliplr(y_coords_quad1(2:end)), y_coords_quad1];
% x_coords_full = [-fliplr(x_coords_quad1(1:end)), x_coords_quad1];
% y_coords_full = [-fliplr(y_coords_quad1(1:end)), y_coords_quad1];%针对起始坐标不是0的情况


% 初始化空的、全尺寸的cell数组
referprojection_full = cell(nx_full, ny_full, nz);
lightMapLibrary_full = cell(nx_full, ny_full, nz);

fprintf('目标维度: %d x %d x %d\n', nx_full, ny_full, nz);

%% 4. 循环遍历旧库，条件性填充到新库
% -------------------------------------------------------------------------
fprintf('正在将第一象限数据对称填充到全尺寸库中 (已修正中心覆盖问题)...\n');

% 计算中心索引
center_idx = nx_quad;

for iz = 1:nz
    for iy = 1:ny_quad
        for ix = 1:nx_quad
            
            % 1. 提取基础数据 (Q1)
            norm_proj_q1 = referprojection_q1{ix, iy, iz};
            lightMap_q1 = lightMapLibrary_q1{ix, iy, iz};
            
            if isempty(norm_proj_q1) || isempty(lightMap_q1)
                continue;
            end
            
            % 2. 准备所有镜像数据 (先计算好，后面按需赋值)
            proj_len_x = size(lightMap_q1, 2); 
            x_proj_q1 = norm_proj_q1(1:proj_len_x);
            y_proj_q1_transposed = norm_proj_q1(proj_len_x+1:end); 
            
            x_proj_flipped = fliplr(x_proj_q1);
            y_proj_flipped_transposed = fliplr(y_proj_q1_transposed);
            
            % Q2 (-x, y)
            lightMap_q2 = fliplr(lightMap_q1); 
            norm_proj_q2 = [x_proj_flipped, y_proj_q1_transposed];
            
            % Q4 (x, -y)
            lightMap_q4 = flipud(lightMap_q1); 
            norm_proj_q4 = [x_proj_q1, y_proj_flipped_transposed];
            
            % Q3 (-x, -y)
            lightMap_q3 = flipud(lightMap_q2); 
            norm_proj_q3 = [x_proj_flipped, y_proj_flipped_transposed];
            
            % --- 计算索引 ---
            % Q1 (+x, +y)
            ix_q1 = center_idx + ix - 1; 
            iy_q1 = center_idx + iy - 1;
            
            % Q2 (-x, +y)
            ix_q2 = center_idx - (ix - 1); 
            iy_q2 = center_idx + iy - 1;
            
            % Q4 (+x, -y)
            ix_q4 = center_idx + ix - 1; 
            iy_q4 = center_idx - (iy - 1);
            
            % Q3 (-x, -y)
            ix_q3 = center_idx - (ix - 1); 
            iy_q3 = center_idx - (iy - 1);
            
            
            % 3. === 核心修改：条件填充逻辑 ===
            
            % [A] 始终填充第一象限 (Q1)
            % 这是最准确的原始数据，永远执行
            referprojection_full{ix_q1, iy_q1, iz} = norm_proj_q1;
            lightMapLibrary_full{ix_q1, iy_q1, iz} = lightMap_q1;
            
            % [B] 仅当不在Y轴上 (x!=0, 即 ix > 1) 时，填充第二象限 (Q2)
            % 这样 ix=1 时就不会用 Q2 覆盖 Q1
            if ix > 1
                referprojection_full{ix_q2, iy_q2, iz} = norm_proj_q2;
                lightMapLibrary_full{ix_q2, iy_q2, iz} = lightMap_q2;
            end
            
            % [C] 仅当不在X轴上 (y!=0, 即 iy > 1) 时，填充第四象限 (Q4)
            % 这样 iy=1 时就不会用 Q4 覆盖 Q1
            if iy > 1
                referprojection_full{ix_q4, iy_q4, iz} = norm_proj_q4;
                lightMapLibrary_full{ix_q4, iy_q4, iz} = lightMap_q4;
            end
            
            % [D] 仅当既不在X轴也不在Y轴 (ix>1 且 iy>1) 时，填充第三象限 (Q3)
            % 避免覆盖轴线上的 Q2 或 Q4
            if ix > 1 && iy > 1
                referprojection_full{ix_q3, iy_q3, iz} = norm_proj_q3;
                lightMapLibrary_full{ix_q3, iy_q3, iz} = lightMap_q3;
            end
            
        end
    end
end
fprintf('数据对称填充完成！\n');

%% 5. 保存结果
% -------------------------------------------------------------------------
referprojection = referprojection_full;
lightMapLibrary = lightMapLibrary_full;
x_coords = x_coords_full;
y_coords = y_coords_full;

save(outputFile, 'referprojection', 'lightMapLibrary', 'x_coords', 'y_coords', 'z_coords');
fprintf('============================================================\n');
fprintf('扩展完成！\n');
fprintf('全尺寸对称参考库已保存到新文件: %s\n', outputFile);




% % % =========================================================================
% % % === 脚本：读取1/4参考库文件，并将其对称扩展为全尺寸库 ===
% % % =========================================================================
% % clear;
% % clc;
% % close all;
% % 
% % %% 1. 定义输入和输出文件
% % % -------------------------------------------------------------------------
% % % 输入文件：您已经生成好的、只包含1/4数据的参考库
% % inputFile = 'MyReferenceLibrary.mat'; 
% % 
% % % 输出文件：扩展后的、完整的、对称填充的参考库
% % outputFile = 'new48.5mmstep1mm.mat';
% % 
% % %% 2. 加载已有的1/4参考库
% % % -------------------------------------------------------------------------
% % fprintf('正在加载已有的1/4参考库: %s ...\n', inputFile);
% % if ~exist(inputFile, 'file')
% %     error('错误：找不到输入文件！请确保 %s 文件与此脚本在同一个目录下。', inputFile);
% % end
% % % 加载后，工作区应包含: referprojection, lightMapLibrary, x_coords, y_coords, z_coords
% % load(inputFile);
% % fprintf('文件加载成功！\n');
% % 
% % %% 3. 根据旧库信息，初始化新的全尺寸库
% % % -------------------------------------------------------------------------
% % % 为了清晰，将从旧库加载的变量重命名
% % referprojection_q1 = referprojection;
% % lightMapLibrary_q1 = lightMapLibrary;
% % x_coords_quad1 = x_coords;
% % y_coords_quad1 = y_coords;
% % % z_coords 保持不变
% % 
% % % 获取第一象限库的维度
% % [nx_quad, ny_quad, nz] = size(referprojection_q1);
% % 
% % % 计算完整库的维度
% % nx_full = nx_quad * 2 - 1;
% % ny_full = ny_quad * 2 - 1;
% % 
% % % 创建代表完整物理范围的新坐标轴
% % x_coords_full = [-fliplr(x_coords_quad1(2:end)), x_coords_quad1];
% % y_coords_full = [-fliplr(y_coords_quad1(2:end)), y_coords_quad1];
% % 
% % % 初始化空的、全尺寸的cell数组
% % referprojection_full = cell(nx_full, ny_full, nz);
% % lightMapLibrary_full = cell(nx_full, ny_full, nz);
% % 
% % fprintf('已初始化全尺寸参考库，目标维度: %d x %d x %d\n', nx_full, ny_full, nz);
% % 
% % %% 4. 循环遍历旧库，将数据对称填充到新库的四个象限
% % % -------------------------------------------------------------------------
% % fprintf('正在将第一象限数据对称填充到全尺寸库中...\n');
% % for iz = 1:nz
% %     for iy = 1:ny_quad
% %         for ix = 1:nx_quad
% % 
% %             % 提取第一象限的基础数据 (Q1)
% %             norm_proj_q1 = referprojection_q1{ix, iy, iz};
% %             lightMap_q1 = lightMapLibrary_q1{ix, iy, iz};
% % 
% %             % 如果当前点为空，则跳过
% %             if isempty(norm_proj_q1) || isempty(lightMap_q1)
% %                 continue;
% %             end
% % 
% %             % --- 关键：进行对称操作，生成其他三个象限的数据 ---
% % 
% %             % 从第一象限的投影中分离出X和Y分量
% %             proj_len_x = size(lightMap_q1, 2); % 探测器X方向像素数
% %             x_proj_q1 = norm_proj_q1(1:proj_len_x);
% %             y_proj_q1_transposed = norm_proj_q1(proj_len_x+1:end); % 这是转置后的Y投影
% % 
% %             % 对投影数据进行镜像
% %             x_proj_flipped = fliplr(x_proj_q1);
% %             y_proj_flipped_transposed = fliplr(y_proj_q1_transposed);
% % 
% %             % 第二象限 (-x, y)
% %             lightMap_q2 = fliplr(lightMap_q1); 
% %             norm_proj_q2 = [x_proj_flipped, y_proj_q1_transposed];
% % 
% %             % 第四象限 (x, -y)
% %             lightMap_q4 = flipud(lightMap_q1); 
% %             norm_proj_q4 = [x_proj_q1, y_proj_flipped_transposed];
% % 
% %             % 第三象限 (-x, -y)
% %             lightMap_q3 = flipud(lightMap_q2); % 在Q2的基础上垂直翻转
% %             norm_proj_q3 = [x_proj_flipped, y_proj_flipped_transposed];
% % 
% %             % --- 计算四个象限在完整库中的索引 ---
% %             center_idx = nx_quad;
% %             ix_q1 = center_idx + ix - 1; iy_q1 = center_idx + iy - 1;
% %             ix_q2 = center_idx - (ix - 1); iy_q2 = center_idx + iy - 1;
% %             ix_q3 = center_idx - (ix - 1); iy_q3 = center_idx - (iy - 1);
% %             ix_q4 = center_idx + ix - 1; iy_q4 = center_idx - (iy - 1);
% % 
% %             % --- 将四个象限的数据全部填充到库中 ---
% %             referprojection_full{ix_q1, iy_q1, iz} = norm_proj_q1; lightMapLibrary_full{ix_q1, iy_q1, iz} = lightMap_q1;
% %             referprojection_full{ix_q2, iy_q2, iz} = norm_proj_q2; lightMapLibrary_full{ix_q2, iy_q2, iz} = lightMap_q2;
% %             referprojection_full{ix_q3, iy_q3, iz} = norm_proj_q3; lightMapLibrary_full{ix_q3, iy_q3, iz} = lightMap_q3;
% %             referprojection_full{ix_q4, iy_q4, iz} = norm_proj_q4; lightMapLibrary_full{ix_q4, iy_q4, iz} = lightMap_q4;
% %         end
% %     end
% % end
% % fprintf('数据对称填充完成！\n');
% % 
% % %% 5. 准备最终保存的变量，并保存到新文件
% % % -------------------------------------------------------------------------
% % % 为了与后续脚本兼容，将最终变量重命名回标准名称
% % referprojection = referprojection_full;
% % lightMapLibrary = lightMapLibrary_full;
% % x_coords = x_coords_full;
% % y_coords = y_coords_full;
% % % z_coords 不需要重命名
% % 
% % save(outputFile, 'referprojection', 'lightMapLibrary', 'x_coords', 'y_coords', 'z_coords');
% % 
% % fprintf('============================================================\n');
% % fprintf('扩展完成！\n');
% % fprintf('全尺寸对称参考库已保存到新文件: %s\n', outputFile);
% 
% % =========================================================================
% % === 脚本：读取1/4参考库文件，并根据特定对称规则扩展为全尺寸库 ===
% % =========================================================================
% clear;
% clc;
% close all;
% 
% %% 1. 定义输入和输出文件
% % -------------------------------------------------------------------------
% inputFile = 'MyReferenceLibrary.mat'; 
% outputFile = 'new48.5mmstep2.5mm_Symmetric.mat'; % 建议修改文件名以免覆盖
% 
% %% 2. 加载已有的1/4参考库
% % -------------------------------------------------------------------------
% fprintf('正在加载已有的1/4参考库: %s ...\n', inputFile);
% if ~exist(inputFile, 'file')
%     error('错误：找不到输入文件！请确保 %s 文件与此脚本在同一个目录下。', inputFile);
% end
% load(inputFile);
% fprintf('文件加载成功！\n');
% 
% %% 3. 根据旧库信息，初始化新的全尺寸库
% % -------------------------------------------------------------------------
% referprojection_q1 = referprojection;
% lightMapLibrary_q1 = lightMapLibrary;
% x_coords_quad1 = x_coords;
% y_coords_quad1 = y_coords;
% 
% % 获取第一象限库的维度
% [nx_quad, ny_quad, nz] = size(referprojection_q1);
% 
% % 计算完整库的维度 (中心点不重复计算，所以是 2*N - 1)
% nx_full = nx_quad * 2 - 1;
% ny_full = ny_quad * 2 - 1;
% 
% % 创建代表完整物理范围的新坐标轴
% x_coords_full = [-fliplr(x_coords_quad1(2:end)), x_coords_quad1];
% y_coords_full = [-fliplr(y_coords_quad1(2:end)), y_coords_quad1];
% 
% % 初始化空的、全尺寸的cell数组
% referprojection_full = cell(nx_full, ny_full, nz);
% lightMapLibrary_full = cell(nx_full, ny_full, nz);
% 
% fprintf('已初始化全尺寸参考库，目标维度: %d x %d x %d\n', nx_full, ny_full, nz);
% 
% %% 4. 循环遍历旧库，根据对称规则填充新库
% % -------------------------------------------------------------------------
% fprintf('正在根据对称规则填充全尺寸库...\n');
% 
% % 定义中心索引位置
% center_idx_x = nx_quad;
% center_idx_y = ny_quad;
% 
% for iz = 1:nz
%     for iy = 1:ny_quad
%         for ix = 1:nx_quad
% 
%             % 提取第一象限的基础数据 (Q1)
%             norm_proj_q1 = referprojection_q1{ix, iy, iz};
%             lightMap_q1 = lightMapLibrary_q1{ix, iy, iz};
% 
%             % 如果当前点为空，则跳过
%             if isempty(norm_proj_q1) || isempty(lightMap_q1)
%                 continue;
%             end
% 
%             % --- 准备基础变换数据 ---
% 
%             % 分离投影数据
%             proj_len_x = size(lightMap_q1, 2); 
%             x_proj_q1 = norm_proj_q1(1:proj_len_x);
%             y_proj_q1_transposed = norm_proj_q1(proj_len_x+1:end); 
% 
%             % 镜像投影数据
%             x_proj_flipped = fliplr(x_proj_q1);
%             y_proj_flipped_transposed = fliplr(y_proj_q1_transposed);
% 
%             % --- 计算全尺寸库中的绝对索引 ---
% 
%             % 第一象限/正半轴索引 (x, y)
%             idx_x_pos = center_idx_x + (ix - 1); 
%             idx_y_pos = center_idx_y + (iy - 1);
% 
%             % 第三象限/负半轴索引 (-x, -y)
%             idx_x_neg = center_idx_x - (ix - 1); 
%             idx_y_neg = center_idx_y - (iy - 1);
% 
%             % --- 逻辑判断与填充 ---
% 
%             if ix == 1 && iy == 1
%                 % =========================================================
%                 % 情况 1: 中心原点 (ix=1 且 iy=1)
%                 % 规则: 只填充第一象限位置 (即中心点本身)
%                 % =========================================================
%                 referprojection_full{idx_x_pos, idx_y_pos, iz} = norm_proj_q1;
%                 lightMapLibrary_full{idx_x_pos, idx_y_pos, iz} = lightMap_q1;
% 
%             elseif ix == 1 || iy == 1
%                 % =========================================================
%                 % 情况 2: 坐标轴上的点 (ix=1 或 iy=1，但不同时为1)
%                 % 规则: 填充第一象限 (正轴) 和 第三象限 (负轴) -> 中心对称
%                 % =========================================================
% 
%                 % 1. 填充第一象限位置 (正轴)
%                 referprojection_full{idx_x_pos, idx_y_pos, iz} = norm_proj_q1;
%                 lightMapLibrary_full{idx_x_pos, idx_y_pos, iz} = lightMap_q1;
% 
%                 % 2. 填充第三象限位置 (负轴，中心对称)
%                 % 数据变换: LightMap 旋转180度 (上下翻转+左右翻转)，投影全部翻转
%                 lightMap_q3 = flipud(fliplr(lightMap_q1)); 
%                 norm_proj_q3 = [x_proj_flipped, y_proj_flipped_transposed];
% 
%                 referprojection_full{idx_x_neg, idx_y_neg, iz} = norm_proj_q3;
%                 lightMapLibrary_full{idx_x_neg, idx_y_neg, iz} = lightMap_q3;
% 
%             else
%                 % =========================================================
%                 % 情况 3: 普通区域点 (ix > 1 且 iy > 1)
%                 % 规则: 填充所有四个象限
%                 % =========================================================
% 
%                 % Q1: 第一象限 (x, y)
%                 referprojection_full{idx_x_pos, idx_y_pos, iz} = norm_proj_q1;
%                 lightMapLibrary_full{idx_x_pos, idx_y_pos, iz} = lightMap_q1;
% 
%                 % Q2: 第二象限 (-x, y) -> 沿Y轴镜像
%                 lightMap_q2 = fliplr(lightMap_q1); 
%                 norm_proj_q2 = [x_proj_flipped, y_proj_q1_transposed];
%                 referprojection_full{idx_x_neg, idx_y_pos, iz} = norm_proj_q2; 
%                 lightMapLibrary_full{idx_x_neg, idx_y_pos, iz} = lightMap_q2;
% 
%                 % Q3: 第三象限 (-x, -y) -> 中心对称
%                 lightMap_q3 = flipud(fliplr(lightMap_q1)); 
%                 norm_proj_q3 = [x_proj_flipped, y_proj_flipped_transposed];
%                 referprojection_full{idx_x_neg, idx_y_neg, iz} = norm_proj_q3; 
%                 lightMapLibrary_full{idx_x_neg, idx_y_neg, iz} = lightMap_q3;
% 
%                 % Q4: 第四象限 (x, -y) -> 沿X轴镜像
%                 lightMap_q4 = flipud(lightMap_q1); 
%                 norm_proj_q4 = [x_proj_q1, y_proj_flipped_transposed];
%                 referprojection_full{idx_x_pos, idx_y_neg, iz} = norm_proj_q4; 
%                 lightMapLibrary_full{idx_x_pos, idx_y_neg, iz} = lightMap_q4;
% 
%             end
%         end
%     end
% end
% fprintf('数据对称填充完成！\n');
% 
% %% 5. 准备最终保存的变量，并保存到新文件
% % -------------------------------------------------------------------------
% referprojection = referprojection_full;
% lightMapLibrary = lightMapLibrary_full;
% x_coords = x_coords_full;
% y_coords = y_coords_full;
% % z_coords 保持不变
% 
% save(outputFile, 'referprojection', 'lightMapLibrary', 'x_coords', 'y_coords', 'z_coords');
% 
% fprintf('============================================================\n');
% fprintf('扩展完成！\n');

% fprintf('全尺寸对称参考库已保存到新文件: %s\n', outputFile);
