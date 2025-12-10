% =========================================================================
% === 自动构建GATE闪烁探测器参考库(含异常文件汇总功能) ===
% =========================================================================
clear;
clc;
close all;

%% 1. 读取JSON参数文件
% -------------------------------------------------------------------------
try
    json_text = fileread('parameters.json');
    params = jsondecode(json_text);
    fprintf('成功读取参数文件: parameters.json\n');
catch ME
    error('无法读取或解析 parameters.json 文件: %s', ME.message);
end

%% 2. 从参数结构体中提取所需变量
% -------------------------------------------------------------------------
% --- 输入/输出设置 ---
inputDirectory = params.io_parameters.output_directory;
outputLibraryFile = params.io_parameters.library_filename;

% --- 探测器参数 ---
det_params = params.detector_parameters;
detector_z_position = det_params.z_position;
detector_pixels_xy = [det_params.pixels_x, det_params.pixels_y];
pixel_size_mm = det_params.pixel_size_mm;

% --- 仿真网格定义 ---
grid_params = params.grid_parameters;
x_coords = grid_params.x.start:grid_params.x.step:grid_params.x.end;
y_coords = grid_params.y.start:grid_params.y.step:grid_params.y.end;
z_coords = grid_params.z.start:grid_params.z.step:grid_params.z.end;

%% 3. 初始化参考库数据结构 
% -------------------------------------------------------------------------
nx = length(x_coords);
ny = length(y_coords);
nz = length(z_coords);

referprojection = cell(nx, ny, nz);
lightMapLibrary = cell(nx, ny, nz); 

% --- !!! 新增：初始化异常文件列表 !!! ---
skipped_files_list = {}; 

fprintf('参考库初始化完成，总网格点数: %d x %d x %d = %d\n', nx, ny, nz, nx*ny*nz);
fprintf('============================================================\n');

%% 4. 循环遍历所有网格点，读取并处理数据
% -------------------------------------------------------------------------
total_points = nx * ny * nz;
processed_count = 0;

for iz = 1:nz
    for iy = 1:ny
        for ix = 1:nx
            processed_count = processed_count + 1;
            
            % --- 获取当前网格点的坐标 ---
            current_x = x_coords(ix);
            current_y = y_coords(iy);
            current_z = z_coords(iz);
            
            % 降低打印频率，每处理100个点打印一次进度，或者打印所有点
            % fprintf('--- 正在处理第 %d/%d 个点: (X=%.2f, Y=%.2f, Z=%.2f) ---\n', ...
            %         processed_count, total_points, current_x, current_y, current_z);
            
            % --- 根据坐标构建输入文件名 ---
            x_fn = strrep(sprintf('%.2f', current_x), '.', 'p');
            y_fn = strrep(sprintf('%.2f', current_y), '.', 'p');
            z_fn = strrep(strrep(sprintf('%.2f', current_z), '.', 'p'), '-', 'm');
            
            filename = sprintf('pos_x%s_y%s_z%sHits.dat', x_fn, y_fn, z_fn);
            full_path = fullfile(inputDirectory, filename);
            
            % --- 检查文件是否存在 ---
            if ~exist(full_path, 'file')
                % warning('文件未找到: %s，跳过此点。', full_path);
                % 也可以选择记录文件不存在的情况
                skipped_files_list{end+1} = sprintf('%s (文件不存在)', filename);
                continue; 
            end
            
            % --- 读取并解析文件 ---
            try
                fid = fopen(full_path, 'r');
                % 根据您的数据格式，这里假设是读取Hits文件
                data = textscan(fid, '%*s %*s %*s %*s %*s %*s %*s %f %*s %*s %f %f %f %*s %*s %*s %*s %*s %*s %s %*s %*s');
                fclose(fid);
            catch ME
                warning('无法读取文件 %s: %s，跳过此点。', full_path, ME.message);
                skipped_files_list{end+1} = sprintf('%s (读取错误)', filename);
                continue;
            end
            
            positions = [data{2}, data{3}, data{4}];
            processNames = data{5};
            
            % --- 筛选光学光子 ---
            tolerance = 0.01;
            detection_indices = strcmp('Transportation', processNames) & (abs(positions(:, 3) - detector_z_position) < tolerance);
            detected_photons_pos = positions(detection_indices, 1:2);
            
            % --- !!! 关键修改：记录没有探测到光子的文件 !!! ---
            if isempty(detected_photons_pos)
                warning('在文件 %s 中没有找到探测到的光子，跳过此点。', filename);
                skipped_files_list{end+1} = filename; % 将文件名加入列表
                continue;
            end
            
            % --- 生成光分布图和投影 ---
            [lightMap, x_projection, y_projection] = createLightMap(detected_photons_pos, detector_pixels_xy, pixel_size_mm);
            
            % --- 存储投影数据 ---
            combined_projection = [x_projection, y_projection']; 
            total_counts = sum(combined_projection);
            
            if total_counts > 0
                normalized_projection = combined_projection / total_counts;
            else
                normalized_projection = zeros(size(combined_projection));
            end
            
            referprojection{ix, iy, iz} = normalized_projection; 
            
            % --- 存储光分布图 ---
            lightMapLibrary{ix, iy, iz} = lightMap;
            
            if mod(processed_count, 100) == 0
                 fprintf('进度: %d/%d (当前探测到 %d 光子)\n', processed_count, total_points, size(detected_photons_pos, 1));
            end
            
        end
    end
end

%% 5. 保存完整的参考库到 .mat 文件
% -------------------------------------------------------------------------
save(outputLibraryFile, 'referprojection','lightMapLibrary', 'x_coords', 'y_coords', 'z_coords');
fprintf('============================================================\n');
fprintf('所有网格点处理完毕！\n');
fprintf('参考库已完整保存到文件: %s\n', outputLibraryFile);

%% 6. !!! 新增：打印异常文件汇总 !!!
% -------------------------------------------------------------------------
fprintf('\n============================================================\n');
fprintf('异常文件汇总 (跳过处理的文件):\n');
if isempty(skipped_files_list)
    fprintf('  无。所有存在的文件均包含有效光子。\n');
else
    fprintf('  共跳过 %d 个文件 (无有效探测光子或文件丢失):\n', length(skipped_files_list));
    fprintf('------------------------------------------------------------\n');
    for k = 1:length(skipped_files_list)
        fprintf('  %d. %s\n', k, skipped_files_list{k});
    end
    fprintf('------------------------------------------------------------\n');
    % 可选：将列表保存到txt文件
    % fid_log = fopen('skipped_files_log.txt', 'w');
    % for k = 1:length(skipped_files_list)
    %     fprintf(fid_log, '%s\n', skipped_files_list{k});
    % end
    % fclose(fid_log);
    % fprintf('  (异常文件列表已保存至 skipped_files_log.txt)\n');
end
fprintf('============================================================\n');


%% 本地函数 (Local Function)
% =========================================================================
function [plane, xcounts, ycounts] = createLightMap(xy_positions, pixel_counts, pixel_dim)
    % 1. 定义物理边界 (Edges)
    % 计算探测器的物理半宽 
    half_width_x = (pixel_counts(1) * pixel_dim) / 2;
    half_width_y = (pixel_counts(2) * pixel_dim) / 2;
    
    % 生成网格边缘向量 (linspace 生成的是边界线，比像素数多1)
    x_edges = linspace(-half_width_x, half_width_x, pixel_counts(1) + 1);
    y_edges = linspace(-half_width_y, half_width_y, pixel_counts(2) + 1);
    
    % 2. 使用二维直方图统计 (使用 histcounts2 进行快速2D分箱)
    % 注意：histcounts2 的 Y 对应行，X 对应列
    h = histcounts2(xy_positions(:,2), xy_positions(:,1), y_edges, x_edges);
    
    plane = uint32(h);
    
    % 3. 计算投影
    xcounts = sum(plane, 1);
    ycounts = sum(plane, 2);
end    