function result = IMAGE_run_recon(cfg)
%IMAGE_RUN_RECON Run IMAGE reconstruction from an app-supplied configuration.
%
% This file is derived from IMAGE_main_recon.m. The original script is kept
% intact; this function form lets the UI build cfg programmatically.

if nargin < 1 || isempty(cfg)
    ui_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(ui_dir);
    cfg = IMAGE_default_config(repo_root);
end

if ~isfield(cfg, 'script_dir') || isempty(cfg.script_dir)
    cfg.script_dir = fileparts(mfilename('fullpath'));
end
method_dir = fullfile(cfg.script_dir, 'method');
if exist(method_dir, 'dir')
    addpath(method_dir);
end

result = struct();
result.recon_images = {};
result.recon_data = {};
result.output_root = '';
result.run_summary = table();

cfg.postprocess_mode = lower(strtrim(char(cfg.postprocess_mode)));
use_selected_files = isfield(cfg, 'selected_files') && ~isempty(cfg.selected_files);
if strcmp(cfg.postprocess_mode, 'single_hole_scan') && ~use_selected_files
    cfg.file_basename = 'calib_events_p%04d.mat';
    cfg.file_indices_to_process = cfg.single_hole_scan.file_indices;
end
cfg.library_file = fullfile(cfg.library_folder, cfg.lse_library_filename);
cfg.render.colorbar_label = get_colorbar_label(cfg.render.colorbar_mode);
cfg.grid = make_reconstruction_grid(cfg);
cfg = configure_method_calibration_paths(cfg);
dirs = prepare_output_dirs(cfg);
result.output_root = dirs.root;

calib = load_calibration_assets(cfg);
locator = prepare_image_locator(cfg);

single_hole_state = [];
if strcmp(cfg.postprocess_mode, 'single_hole')
    single_hole_state = init_single_hole_state(cfg);
end
single_hole_scan_state = [];
if strcmp(cfg.postprocess_mode, 'single_hole_scan')
    single_hole_scan_state = init_single_hole_scan_state(cfg);
end

run_summary = cell(0, 7);

if use_selected_files
    files_to_process = cfg.selected_files(:)';
else
    files_to_process = cell(1, numel(cfg.file_indices_to_process));
    for kk = 1:numel(cfg.file_indices_to_process)
        files_to_process{kk} = fullfile(cfg.input_data_folder, ...
            sprintf(cfg.file_basename, cfg.file_indices_to_process(kk)));
    end
end

for k = 1:numel(files_to_process)
    input_file = files_to_process{k};
    [~, current_stem, current_ext] = fileparts(input_file);
    current_filename = [current_stem, current_ext];
    if use_selected_files
        file_idx = infer_file_index(current_filename, k);
        file_tag = sanitize_file_tag(current_stem);
    else
        file_idx = cfg.file_indices_to_process(k);
        file_tag = sprintf('p%04d', file_idx);
    end
    if ~exist(input_file, 'file')
        fprintf('Skip missing file: %s\n', input_file);
        continue;
    end

    fprintf('\nProcessing %s [%s]...\n', current_filename, cfg.postprocess_mode);
    t_file_start = tic;

    data_struct = load(input_file);
    if ~isfield(data_struct, 'planeset')
        warning('File does not contain planeset: %s', input_file);
        continue;
    end

    planeset = data_struct.planeset;
    if isempty(planeset)
        fprintf('Skip empty planeset: %s\n', current_filename);
        continue;
    end

    localized = locate_image_events(planeset, locator, cfg);
    localized.file_idx = file_idx;
    localized.filename = current_filename;
    localized.algorithm_name = locator.algorithm_name;

    localized = apply_quality_filter(localized, cfg);
    localized = apply_spatial_correction(localized, calib, cfg);
    localized.event_weights = compute_event_weights(localized.true_x_mm, localized.true_y_mm, calib, cfg);

    recon = reconstruct_image(localized.true_x_mm, localized.true_y_mm, calib, cfg);

    recon_tag = sprintf('%s_%s_%s', locator.algorithm_name, cfg.postprocess_mode, file_tag);
    recon_save_path = fullfile(dirs.reconstruction, sprintf('Recon_%s.png', recon_tag));
    plot_reconstruction_image(recon.display_counts, cfg.grid.bin_centers, cfg.grid.axis_limit, ...
        recon_save_path, cfg.render);

    recon_data_path = fullfile(dirs.reconstruction, sprintf('ReconData_%s.mat', recon_tag));
    save(recon_data_path, 'localized', 'recon', 'cfg', '-v7.3');
    result.recon_images{end+1} = recon_save_path; %#ok<AGROW>
    result.recon_data{end+1} = recon_data_path; %#ok<AGROW>
    fprintf('Saved reconstruction image: %s\n', recon_save_path);
    fprintf('Saved reconstruction data:  %s\n', recon_data_path);

    switch cfg.postprocess_mode
        case 'flood'
            postprocess_flood(recon, cfg, dirs, file_idx, file_tag, locator.algorithm_name);
        case 'gridmask'
            postprocess_gridmask(recon, cfg, dirs, file_idx, file_tag, locator.algorithm_name);
        case 'slit'
            postprocess_slit(localized, recon, cfg, dirs, file_idx, file_tag, locator.algorithm_name);
        case 'single_hole'
            single_hole_state = accumulate_single_hole(single_hole_state, recon, cfg, dirs, file_idx, file_tag, locator.algorithm_name);
        case 'single_hole_scan'
            single_hole_scan_state = accumulate_single_hole_scan(single_hole_scan_state, ...
                localized, recon, cfg, dirs, file_idx, file_tag, locator.algorithm_name);
        otherwise
            error('Unknown postprocess_mode: %s', cfg.postprocess_mode);
    end

    run_summary(end+1, :) = {current_filename, file_idx, file_tag, numel(planeset), ...
        numel(localized.true_x_mm), cfg.postprocess_mode, toc(t_file_start)}; %#ok<SAGROW>
    fprintf('Finished %s in %.2f s.\n', current_filename, toc(t_file_start));
end

if strcmp(cfg.postprocess_mode, 'single_hole')
    finalize_single_hole(single_hole_state, cfg, dirs, locator.algorithm_name);
end
if strcmp(cfg.postprocess_mode, 'single_hole_scan')
    finalize_single_hole_scan(single_hole_scan_state, cfg, dirs, locator.algorithm_name);
end

if ~isempty(run_summary)
    RunSummary = cell2table(run_summary, 'VariableNames', ...
        {'Filename', 'FileIndex', 'FileTag', 'NumEvents', 'NumValidEvents', 'PostprocessMode', 'ElapsedSeconds'});
    writetable(RunSummary, fullfile(dirs.root, sprintf('RunSummary_%s.csv', cfg.postprocess_mode)));
    save(fullfile(dirs.root, sprintf('RunSummary_%s.mat', cfg.postprocess_mode)), 'RunSummary');
    result.run_summary = RunSummary;
end

fprintf('\nIMAGE pipeline complete.\n');
end

%% ========================================================================
% Configuration and setup
% =========================================================================
function cfg = default_image_config(script_dir, repo_root)
    cfg.script_dir = script_dir;
    cfg.repo_root = repo_root;
    cfg.input_data_folder = fullfile(repo_root, 'input_data');
    cfg.output_root = fullfile(script_dir, 'IMAGE_Recon_Results');

    cfg.calib_root = fullfile(repo_root, 'CALIB');
    cfg.lse_calibration_folder = fullfile(cfg.calib_root, ...
        'Calibration_Results_FloodCDF_LSE');
    cfg.anger_calibration_folder = fullfile(cfg.calib_root, ...
        'Calibration_Results_AngerCDF');
    cfg.spatial_cdf_path = fullfile(cfg.lse_calibration_folder, 'CDF_Spatial_LUT.mat');
    cfg.ucm_path = fullfile(cfg.lse_calibration_folder, 'Uniformity_UCM.mat');
    cfg.library_folder = fullfile(repo_root, 'library-construction');
    cfg.lse_library_filename = 'newSizeDetgamma0.2mmfull.mat';
    cfg.library_file = fullfile(cfg.library_folder, cfg.lse_library_filename);

    cfg.lib_phys_start = -25.0;
    cfg.lib_phys_end = 25.0;
    cfg.show_size_mm = 50.0;
    cfg.display_resolution = 512;
    cfg.detector_pitch_mm = 6.0;
end

function cfg = configure_method_calibration_paths(cfg)
    method = lower(strtrim(char(cfg.localization.method)));
    switch method
        case {'anger', 'anger_standard', 'anger_rtp'}
            cfg.spatial_cdf_path = fullfile(cfg.anger_calibration_folder, 'CDF_Spatial_LUT.mat');
            cfg.ucm_path = fullfile(cfg.anger_calibration_folder, 'Uniformity_UCM.mat');
            cfg.calibration_source = 'AngerCDF';
        otherwise
            cfg.spatial_cdf_path = fullfile(cfg.lse_calibration_folder, 'CDF_Spatial_LUT.mat');
            cfg.ucm_path = fullfile(cfg.lse_calibration_folder, 'Uniformity_UCM.mat');
            cfg.calibration_source = 'FloodCDF_LSE';
    end
end

function grid = make_reconstruction_grid(cfg)
    grid.axis_limit = cfg.show_size_mm / 2;
    grid.bin_edges = linspace(-grid.axis_limit, grid.axis_limit, cfg.display_resolution + 1);
    grid.bin_centers = (grid.bin_edges(1:end-1) + grid.bin_edges(2:end)) / 2;
    grid.pixel_area_mm2 = (cfg.show_size_mm / cfg.display_resolution)^2;
    grid.bin_step_mm = mean(diff(grid.bin_centers));
end

function dirs = prepare_output_dirs(cfg)
    dirs.root = cfg.output_root;
    dirs.reconstruction = fullfile(dirs.root, 'reconstruction');
    dirs.flood = fullfile(dirs.root, 'flood');
    dirs.gridmask = fullfile(dirs.root, 'gridmask');
    dirs.slit = fullfile(dirs.root, 'slit');
    dirs.single_hole = fullfile(dirs.root, 'single_hole');
    dirs.single_hole_scan = fullfile(dirs.root, 'single_hole_scan');
    dirs.single_hole_srfigure = fullfile(dirs.single_hole_scan, 'SRfigure');

    ensure_dir(dirs.root);
    ensure_dir(dirs.reconstruction);
    ensure_dir(dirs.flood);
    ensure_dir(fullfile(dirs.flood, 'projections'));
    ensure_dir(dirs.gridmask);
    ensure_dir(fullfile(dirs.gridmask, 'profiles'));
    ensure_dir(fullfile(dirs.gridmask, 'tables'));
    ensure_dir(dirs.slit);
    ensure_dir(fullfile(dirs.slit, 'aligned_profiles'));
    ensure_dir(dirs.single_hole);
    ensure_dir(fullfile(dirs.single_hole, 'per_file'));
    ensure_dir(fullfile(dirs.single_hole, 'composite'));
    ensure_dir(dirs.single_hole_scan);
    ensure_dir(fullfile(dirs.single_hole_scan, 'per_file'));
    ensure_dir(fullfile(dirs.single_hole_scan, 'figures'));
    ensure_dir(fullfile(dirs.single_hole_scan, 'tables'));
    ensure_dir(dirs.single_hole_srfigure);
end

function ensure_dir(path_name)
    if ~exist(path_name, 'dir')
        mkdir(path_name);
    end
end

function file_idx = infer_file_index(filename, fallback_idx)
    token = regexp(filename, 'p(?:gridmask|slitx|slity|flood)?(\d+)', 'tokens', 'once');
    if isempty(token)
        file_idx = fallback_idx;
        return;
    end
    file_idx = str2double(token{1});
    if ~isfinite(file_idx)
        file_idx = fallback_idx;
    end
end

function tag = sanitize_file_tag(raw_tag)
    tag = regexprep(char(raw_tag), '[^\w\-.]', '_');
    if isempty(tag)
        tag = 'selected_file';
    end
end

function tag = summarize_file_tags(file_tags)
    file_tags = file_tags(:)';
    if numel(file_tags) == 1
        tag = file_tags{1};
    elseif numel(file_tags) == 2
        tag = sprintf('%s_%s', file_tags{1}, file_tags{2});
    else
        tag = sprintf('%s_to_%s_N%d', file_tags{1}, file_tags{end}, numel(file_tags));
    end
    tag = sanitize_file_tag(tag);
end

function calib = load_calibration_assets(cfg)
    calib = struct();
    calib.has_cdf = false;
    calib.has_ucm = false;

    if cfg.enable_ucm_correction && ~cfg.enable_cdf_correction
        error('UCM correction depends on the CDF-corrected coordinate grid.');
    end

    if cfg.enable_cdf_correction
        if ~exist(cfg.spatial_cdf_path, 'file')
            error('CDF LUT not found: %s', cfg.spatial_cdf_path);
        end
        data = load(cfg.spatial_cdf_path, 'cdf_centers', 'mapped_X_phys', 'mapped_Y_phys');
        calib.cdf_centers = data.cdf_centers;
        calib.mapped_X_phys = data.mapped_X_phys;
        calib.mapped_Y_phys = data.mapped_Y_phys;
        calib.has_cdf = true;
    end

    if cfg.enable_ucm_correction
        if ~exist(cfg.ucm_path, 'file')
            error('UCM file not found: %s', cfg.ucm_path);
        end
        data = load(cfg.ucm_path, 'UCM', 'ufov_mask');
        calib.UCM = data.UCM;
        if isfield(data, 'ufov_mask')
            calib.ufov_mask = data.ufov_mask;
        else
            calib.ufov_mask = true(size(data.UCM));
        end
        calib.has_ucm = true;
    end
end

%% ========================================================================
% Position filtering and calibration
% =========================================================================
function localized = apply_quality_filter(localized, cfg)
    if isfield(localized, 'valid_mask') && numel(localized.valid_mask) == numel(localized.raw_x_mm)
        base_valid = localized.valid_mask(:);
    else
        base_valid = true(size(localized.raw_x_mm));
    end
    base_valid = base_valid & isfinite(localized.raw_x_mm) & isfinite(localized.raw_y_mm);

    if cfg.localization.enable_quality_filter && any(base_valid & isfinite(localized.residual))
        finite_residual = localized.residual(base_valid & isfinite(localized.residual));
        threshold = prctile(finite_residual, cfg.localization.residual_percentile_cutoff);
        localized.valid_mask = base_valid & isfinite(localized.residual) & localized.residual <= threshold;
    else
        localized.valid_mask = base_valid;
    end

    localized.raw_x_mm = localized.raw_x_mm(localized.valid_mask);
    localized.raw_y_mm = localized.raw_y_mm(localized.valid_mask);
    localized.residual = localized.residual(localized.valid_mask);
end

function localized = apply_spatial_correction(localized, calib, cfg)
    if cfg.enable_cdf_correction && calib.has_cdf
        localized.true_x_mm = interp1(calib.cdf_centers, calib.mapped_X_phys, ...
            double(localized.raw_x_mm), 'linear', 'extrap');
        localized.true_y_mm = interp1(calib.cdf_centers, calib.mapped_Y_phys, ...
            double(localized.raw_y_mm), 'linear', 'extrap');
    else
        localized.true_x_mm = double(localized.raw_x_mm);
        localized.true_y_mm = double(localized.raw_y_mm);
    end
end

function event_weights = compute_event_weights(true_x, true_y, calib, cfg)
    if cfg.enable_ucm_correction && calib.has_ucm
        idx_x = round((true_x - cfg.grid.bin_centers(1)) / cfg.grid.bin_step_mm) + 1;
        idx_y = round((true_y - cfg.grid.bin_centers(1)) / cfg.grid.bin_step_mm) + 1;
        idx_x = max(1, min(cfg.display_resolution, idx_x));
        idx_y = max(1, min(cfg.display_resolution, idx_y));
        event_weights = calib.UCM(sub2ind(size(calib.UCM), idx_y, idx_x));
    else
        event_weights = ones(size(true_x));
    end
    event_weights = double(event_weights(:));
end

%% ========================================================================
% Reconstruction and plotting
% =========================================================================
function recon = reconstruct_image(true_x, true_y, calib, cfg)
    raw_counts = histcounts2(true_y, true_x, cfg.grid.bin_edges, cfg.grid.bin_edges);

    if cfg.render.enable_gaussian_diffusion
        smoothed_counts = gaussian_smooth_counts(raw_counts, ...
            cfg.render.gaussian_sigma_bins, cfg.render.gaussian_radius_bins, ...
            cfg.render.gaussian_kernel_normalization);
    else
        smoothed_counts = raw_counts;
    end

    if cfg.enable_ucm_correction && calib.has_ucm
        analysis_counts = smoothed_counts .* calib.UCM;
        if isfield(calib, 'ufov_mask')
            analysis_counts(~calib.ufov_mask) = 0;
        end
    else
        analysis_counts = smoothed_counts;
    end

    display_counts = convert_colorbar_units(analysis_counts, cfg.grid.pixel_area_mm2, cfg.render.colorbar_mode);

    recon.raw_counts = raw_counts;
    recon.smoothed_counts = smoothed_counts;
    recon.analysis_counts = analysis_counts;
    recon.display_counts = display_counts;
end

function smoothed_counts = gaussian_smooth_counts(counts, sigma_bins, radius_bins, normalization_mode)
    if sigma_bins <= 0
        smoothed_counts = counts;
        return;
    end

    radius_bins = max(1, round(radius_bins));
    [kx, ky] = meshgrid(-radius_bins:radius_bins, -radius_bins:radius_bins);
    kernel = exp(-(kx.^2 + ky.^2) ./ (2 * sigma_bins^2));

    switch lower(normalization_mode)
        case 'sum1'
            kernel = kernel ./ (sum(kernel(:)) + eps);
        case 'peak1'
            kernel = kernel ./ (max(kernel(:)) + eps);
        otherwise
            error('Unknown gaussian_kernel_normalization: %s', normalization_mode);
    end

    smoothed_counts = conv2(double(counts), kernel, 'same');
end

function display_counts = convert_colorbar_units(counts, pixel_area_mm2, colorbar_mode)
    switch lower(colorbar_mode)
        case 'counts'
            display_counts = counts;
        case 'counts_per_mm2'
            display_counts = counts ./ pixel_area_mm2;
        otherwise
            error('Unknown colorbar_mode: %s', colorbar_mode);
    end
end

function label = get_colorbar_label(colorbar_mode)
    switch lower(colorbar_mode)
        case 'counts'
            label = 'Counts / pixel';
        case 'counts_per_mm2'
            label = 'Counts / mm^2';
        otherwise
            error('Unknown colorbar_mode: %s', colorbar_mode);
    end
end

function plot_reconstruction_image(img_data, bin_centers, axis_limit, save_path, render_cfg)
    fig = figure('Visible', 'off', 'Position', [100, 100, 650, 600], 'Color', 'w');
    imagesc(bin_centers, bin_centers, img_data);
    set(gca, 'YDir', 'normal');
    colormap(jet(256));
    axis equal tight;
    cb = colorbar;
    ylabel(cb, render_cfg.colorbar_label, 'Interpreter', 'none');
    apply_image_clim(img_data, render_cfg);
    xlabel('X (mm)');
    ylabel('Y (mm)');

    if render_cfg.enable_grid_lines
        grid_ticks = -axis_limit : 5 : axis_limit;
        set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks, ...
            'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
        grid on;
    else
        set(gca, 'XTick', -25:5:25, 'YTick', -25:5:25);
    end
    set(gca, 'Layer', 'top');

    exportgraphics(fig, save_path, 'Resolution', 300);
    close(fig);
end

function apply_image_clim(img_data, render_cfg)
    if ~render_cfg.enable_slit_style_clim
        return;
    end

    pixels = img_data(isfinite(img_data));
    if isempty(pixels)
        min_disp = 0;
        max_disp = 1;
    else
        if render_cfg.clim_low_percentile <= 0
            min_disp = min(pixels(:));
        else
            min_disp = prctile(pixels(:), render_cfg.clim_low_percentile);
        end

        if render_cfg.clim_high_percentile >= 100
            max_disp = max(pixels(:));
        else
            max_disp = prctile(pixels(:), render_cfg.clim_high_percentile);
        end
    end

    if max_disp <= min_disp || isnan(min_disp) || isnan(max_disp)
        min_disp = 0;
        max_disp = max(1, max(img_data(:)));
    end
    clim([min_disp, max_disp]);
end

%% ========================================================================
% Flood post-processing
% =========================================================================
function postprocess_flood(recon, cfg, dirs, ~, file_tag, algo_name)
    flood_dir = fullfile(dirs.flood, 'projections');
    img = recon.analysis_counts;
    x_projection = sum(img, 1);
    y_projection = sum(img, 2)';

    fig = figure('Visible', 'off', 'Position', [100, 100, 1000, 420], 'Color', 'w');
    subplot(1, 2, 1);
    plot(cfg.grid.bin_centers, x_projection, 'b-', 'LineWidth', 1.5);
    xlabel('X (mm)');
    ylabel('Counts');
    grid on;
    title('X projection');

    subplot(1, 2, 2);
    plot(cfg.grid.bin_centers, y_projection, 'r-', 'LineWidth', 1.5);
    xlabel('Y (mm)');
    ylabel('Counts');
    grid on;
    title('Y projection');

    save_path = fullfile(flood_dir, sprintf('FloodProjection_%s_%s.png', algo_name, file_tag));
    exportgraphics(fig, save_path, 'Resolution', 300);
    close(fig);

    ProjectionTable = table(cfg.grid.bin_centers(:), x_projection(:), y_projection(:), ...
        'VariableNames', {'Position_mm', 'XProjection_Counts', 'YProjection_Counts'});
    writetable(ProjectionTable, fullfile(flood_dir, sprintf('FloodProjection_%s_%s.csv', algo_name, file_tag)));
    save(fullfile(flood_dir, sprintf('FloodProjection_%s_%s.mat', algo_name, file_tag)), ...
        'ProjectionTable', 'x_projection', 'y_projection');
end

%% ========================================================================
% Grid-mask post-processing
% =========================================================================
function postprocess_gridmask(recon, cfg, dirs, file_idx, file_tag, algo_name)
    profile_dir = fullfile(dirs.gridmask, 'profiles');
    table_dir = fullfile(dirs.gridmask, 'tables');

    profile = extract_manual_cross_profile(recon.display_counts, cfg, profile_dir, file_idx, file_tag, algo_name);
    if isempty(profile)
        warning('No gridmask profile was selected for p%04d.', file_idx);
        return;
    end

    x_rows = profile_feature_rows('X', cfg.grid.bin_centers, profile.raw_x, profile.smooth_x, cfg.gridmask);
    y_rows = profile_feature_rows('Y', cfg.grid.bin_centers, profile.raw_y, profile.smooth_y, cfg.gridmask);
    rows = [x_rows; y_rows];

    if isempty(rows)
        FeatureTable = cell2table(cell(0, 10), 'VariableNames', gridmask_feature_column_names());
    else
        FeatureTable = cell2table(rows, 'VariableNames', gridmask_feature_column_names());
    end

    csv_path = fullfile(table_dir, sprintf('GridmaskPeakValley_%s_%s.csv', algo_name, file_tag));
    mat_path = fullfile(table_dir, sprintf('GridmaskPeakValley_%s_%s.mat', algo_name, file_tag));
    writetable(FeatureTable, csv_path);
    save(mat_path, 'FeatureTable', 'profile');

    [SpacingTable, SpacingSummary] = compute_gridmask_profile_spacing(FeatureTable, cfg.gridmask);
    spacing_tag = sprintf('%s_%s', algo_name, file_tag);
    writetable(SpacingTable, fullfile(table_dir, sprintf('GridmaskSpacing_%s.csv', spacing_tag)));
    writetable(SpacingSummary, fullfile(table_dir, sprintf('GridmaskSpacingSummary_%s.csv', spacing_tag)));
    save(fullfile(table_dir, sprintf('GridmaskSpacing_%s.mat', spacing_tag)), ...
        'SpacingTable', 'SpacingSummary', 'FeatureTable', 'profile');

    plot_gridmask_profiles(profile, FeatureTable, cfg, profile_dir, file_tag, algo_name);
end

function names = gridmask_feature_column_names()
    names = {'Axis', 'FeatureType', 'FeatureIndex', 'Position_mm', 'Value', ...
        'Sigma_mm', 'FWHM_mm', 'LeftPeakIndex', 'RightPeakIndex', 'Source'};
end

function [SpacingTable, SpacingSummary] = compute_gridmask_profile_spacing(feature_table, gm_cfg)
    spacing_names = {'Axis', 'PairIndex', 'LeftPeakIndex', 'RightPeakIndex', ...
        'LeftPeakPosition_mm', 'RightPeakPosition_mm', 'PairMidpoint_mm', ...
        'MeasuredSpacing_mm', 'NominalSpacing_mm', 'SignedError_mm', ...
        'AbsoluteError_mm', 'EstimatedNominalIntervals', 'EstimatedSkippedPeaks', ...
        'UsedForError'};
    spacing_types = {'cell', 'double', 'double', 'double', 'double', 'double', ...
        'double', 'double', 'double', 'double', 'double', 'double', 'double', 'logical'};
    SpacingTable = table('Size', [0, numel(spacing_names)], ...
        'VariableTypes', spacing_types, 'VariableNames', spacing_names);

    nominal_spacing = gm_cfg.nominal_spacing_mm;
    max_single_spacing = gm_cfg.max_adjacent_spacing_factor * nominal_spacing;
    axes_to_measure = {'X', 'Y'};

    for aa = 1:numel(axes_to_measure)
        axis_name = axes_to_measure{aa};
        peak_mask = strcmp(feature_table.Axis, axis_name) & strcmp(feature_table.FeatureType, 'peak');
        peaks = feature_table(peak_mask, :);
        if height(peaks) < 2
            continue;
        end

        [~, order] = sort(peaks.Position_mm);
        peaks = peaks(order, :);
        pair_count = height(peaks) - 1;
        measured_spacing = diff(peaks.Position_mm);
        nominal_intervals = max(1, round(measured_spacing ./ nominal_spacing));
        skipped_peaks = max(0, nominal_intervals - 1);
        used_for_error = nominal_intervals == 1 & measured_spacing <= max_single_spacing;
        signed_error = measured_spacing - nominal_spacing;
        absolute_error = abs(signed_error);

        rows = table(repmat({axis_name}, pair_count, 1), (1:pair_count)', ...
            peaks.FeatureIndex(1:end-1), peaks.FeatureIndex(2:end), ...
            peaks.Position_mm(1:end-1), peaks.Position_mm(2:end), ...
            (peaks.Position_mm(1:end-1) + peaks.Position_mm(2:end)) ./ 2, ...
            measured_spacing, repmat(nominal_spacing, pair_count, 1), ...
            signed_error, absolute_error, nominal_intervals, skipped_peaks, used_for_error, ...
            'VariableNames', spacing_names);
        SpacingTable = [SpacingTable; rows]; %#ok<AGROW>
    end

    summary_names = {'Axis', 'NumFittedPeaks', 'NumConsecutiveIntervals', ...
        'NumUsedAdjacentPairs', 'NumExcludedIntervals', 'MeanSpacing_mm', ...
        'MeanSignedError_mm', 'MeanAbsoluteError_mm', 'RMSError_mm', ...
        'MaxAbsoluteError_mm'};
    summary_types = {'cell', 'double', 'double', 'double', 'double', ...
        'double', 'double', 'double', 'double', 'double'};
    SpacingSummary = table('Size', [0, numel(summary_names)], ...
        'VariableTypes', summary_types, 'VariableNames', summary_names);

    for aa = 1:numel(axes_to_measure)
        axis_name = axes_to_measure{aa};
        num_peaks = nnz(strcmp(feature_table.Axis, axis_name) & strcmp(feature_table.FeatureType, 'peak'));
        axis_pairs = SpacingTable(strcmp(SpacingTable.Axis, axis_name), :);
        valid_pairs = axis_pairs(axis_pairs.UsedForError, :);
        SpacingSummary = [SpacingSummary; make_spacing_summary_row(axis_name, num_peaks, axis_pairs, valid_pairs)]; %#ok<AGROW>
    end

    all_peaks = nnz(strcmp(feature_table.FeatureType, 'peak'));
    valid_pairs = SpacingTable(SpacingTable.UsedForError, :);
    SpacingSummary = [SpacingSummary; make_spacing_summary_row('Combined', all_peaks, SpacingTable, valid_pairs)];
end

function row = make_spacing_summary_row(axis_name, num_peaks, all_pairs, valid_pairs)
    if isempty(valid_pairs)
        mean_spacing = NaN;
        mean_signed_error = NaN;
        mean_absolute_error = NaN;
        rms_error = NaN;
        max_absolute_error = NaN;
    else
        mean_spacing = mean(valid_pairs.MeasuredSpacing_mm);
        mean_signed_error = mean(valid_pairs.SignedError_mm);
        mean_absolute_error = mean(valid_pairs.AbsoluteError_mm);
        rms_error = sqrt(mean(valid_pairs.SignedError_mm .^ 2));
        max_absolute_error = max(valid_pairs.AbsoluteError_mm);
    end
    row = table({axis_name}, num_peaks, height(all_pairs), height(valid_pairs), ...
        height(all_pairs) - height(valid_pairs), mean_spacing, mean_signed_error, ...
        mean_absolute_error, rms_error, max_absolute_error, 'VariableNames', ...
        {'Axis', 'NumFittedPeaks', 'NumConsecutiveIntervals', ...
        'NumUsedAdjacentPairs', 'NumExcludedIntervals', 'MeanSpacing_mm', ...
        'MeanSignedError_mm', 'MeanAbsoluteError_mm', 'RMSError_mm', ...
        'MaxAbsoluteError_mm'});
end

function profile = extract_manual_cross_profile(img_data, cfg, profile_dir, ~, file_tag, algo_name)
    profile = [];
    fig = figure('Visible', 'on', 'Position', [100, 100, 650, 600], 'Color', 'w');
    imagesc(cfg.grid.bin_centers, cfg.grid.bin_centers, img_data);
    set(gca, 'YDir', 'normal');
    colormap(jet(256));
    axis equal tight;
    cb = colorbar;
    ylabel(cb, cfg.render.colorbar_label, 'Interpreter', 'none');
    apply_image_clim(img_data, cfg.render);
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title({'Gridmask manual profile', 'Click near the cross-profile center'}, 'Interpreter', 'none');
    hold on;

    if cfg.render.enable_grid_lines
        grid_ticks = -cfg.grid.axis_limit : 5 : cfg.grid.axis_limit;
        set(gca, 'XTick', grid_ticks, 'YTick', grid_ticks, ...
            'GridColor', 'w', 'GridAlpha', 0.5, 'GridLineStyle', '-');
        grid on;
    else
        set(gca, 'XTick', -25:5:25, 'YTick', -25:5:25);
    end

    [cx, cy] = ginput(1);
    if isempty(cx) || isempty(cy)
        close(fig);
        return;
    end

    [~, idx_x_click] = min(abs(cfg.grid.bin_centers - cx));
    [~, idx_y_click] = min(abs(cfg.grid.bin_centers - cy));
    search_radius_bins = max(1, round(cfg.gridmask.manual_peak_search_radius_mm / cfg.grid.bin_step_mm));
    x1 = max(1, idx_x_click - search_radius_bins);
    x2 = min(numel(cfg.grid.bin_centers), idx_x_click + search_radius_bins);
    y1 = max(1, idx_y_click - search_radius_bins);
    y2 = min(numel(cfg.grid.bin_centers), idx_y_click + search_radius_bins);

    local_img = img_data(y1:y2, x1:x2);
    [~, local_max_idx] = max(local_img(:));
    [local_y, local_x] = ind2sub(size(local_img), local_max_idx);
    idx_x = x1 + local_x - 1;
    idx_y = y1 + local_y - 1;
    peak_x = cfg.grid.bin_centers(idx_x);
    peak_y = cfg.grid.bin_centers(idx_y);

    plot(cx, cy, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
    plot(cx, cy, 'ko', 'MarkerSize', 15, 'LineWidth', 2);
    xline(peak_x, 'w-', 'LineWidth', 1.2);
    yline(peak_y, 'w-', 'LineWidth', 1.2);
    plot(peak_x, peak_y, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    plot(peak_x, peak_y, 'wo', 'MarkerSize', 15, 'LineWidth', 2);
    drawnow;

    exportgraphics(fig, fullfile(profile_dir, ...
        sprintf('GridmaskProfileLocation_%s_%s.png', algo_name, file_tag)), 'Resolution', 300);
    close(fig);

    profile.raw_x = img_data(idx_y, :);
    profile.raw_y = img_data(:, idx_x)';
    profile.smooth_x = smoothdata(profile.raw_x, 'gaussian', cfg.gridmask.profile_smooth_window_bins);
    profile.smooth_y = smoothdata(profile.raw_y, 'gaussian', cfg.gridmask.profile_smooth_window_bins);
    profile.click_x = cfg.grid.bin_centers(idx_x_click);
    profile.click_y = cfg.grid.bin_centers(idx_y_click);
    profile.peak_x = peak_x;
    profile.peak_y = peak_y;
    profile.idx_x = idx_x;
    profile.idx_y = idx_y;
    profile.value_label = cfg.render.colorbar_label;
end

function rows = profile_feature_rows(axis_name, x, raw_profile, smooth_profile, gm_cfg)
    rows = cell(0, 10);
    y = double(smooth_profile(:))';
    raw_profile = double(raw_profile(:))';

    if isempty(y) || all(~isfinite(y)) || max(y) <= min(y)
        return;
    end

    y_range = max(y) - min(y);
    min_height = min(y) + gm_cfg.min_peak_relative_height * y_range;
    candidate_idx = find(y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
    candidate_idx = candidate_idx(y(candidate_idx) >= min_height);

    if isempty(candidate_idx)
        [~, max_idx] = max(y);
        candidate_idx = max_idx;
    end

    [~, order] = sort(y(candidate_idx), 'descend');
    selected = [];
    for oi = 1:numel(order)
        idx = candidate_idx(order(oi));
        if isempty(selected) || all(abs(x(idx) - x(selected)) >= gm_cfg.min_peak_distance_mm)
            selected(end+1) = idx; %#ok<AGROW>
        end
        if numel(selected) >= gm_cfg.max_num_peaks
            break;
        end
    end

    selected = sort(selected);
    peak_positions = zeros(numel(selected), 1);

    for pp = 1:numel(selected)
        fit = fit_local_gaussian_peak(x, raw_profile, selected(pp), gm_cfg);
        peak_positions(pp) = fit.center_mm;
        rows(end+1, :) = {axis_name, 'peak', pp, fit.center_mm, fit.peak_value, ...
            fit.sigma_mm, fit.fwhm_mm, NaN, NaN, 'gaussian_fit'}; %#ok<AGROW>
    end

    for vv = 1:max(0, numel(peak_positions)-1)
        left_peak = peak_positions(vv);
        right_peak = peak_positions(vv + 1);
        valley_mask = x >= left_peak & x <= right_peak;
        if any(valley_mask)
            valley_x = x(valley_mask);
            valley_y = y(valley_mask);
            [valley_value, min_idx] = min(valley_y);
            valley_pos = valley_x(min_idx);
            rows(end+1, :) = {axis_name, 'valley', vv, valley_pos, valley_value, ...
                NaN, NaN, vv, vv + 1, 'between_fitted_peaks'}; %#ok<AGROW>
        end
    end
end

function fit = fit_local_gaussian_peak(x, y, peak_idx, gm_cfg)
    center0 = x(peak_idx);
    fit_mask = abs(x - center0) <= gm_cfg.fit_half_width_mm;
    if nnz(fit_mask) < 5
        fit_mask = abs((1:numel(x)) - peak_idx) <= 3;
    end

    xw = double(x(fit_mask));
    yw = double(y(fit_mask));
    bg0 = min(yw);
    amp0 = max(yw) - bg0;
    if amp0 <= 0
        amp0 = max(yw);
    end
    p0 = [amp0, center0, gm_cfg.initial_sigma_mm, bg0];

    obj = @(p) gaussian_sse_with_penalty(p, xw, yw, min(xw), max(xw));
    options = optimset('Display', 'off');
    p_opt = fminsearch(obj, p0, options);

    fit.amp = p_opt(1);
    fit.center_mm = p_opt(2);
    fit.sigma_mm = abs(p_opt(3));
    fit.background = p_opt(4);
    fit.peak_value = p_opt(1) + p_opt(4);
    fit.fwhm_mm = 2.3548 * fit.sigma_mm;
end

function err = gaussian_sse_with_penalty(p, x, y, xmin, xmax)
    sigma = abs(p(3));
    model = gaussian_1d(p, x);
    err = sum((y - model).^2);
    if sigma <= eps || p(2) < xmin || p(2) > xmax || p(1) < 0
        err = err + 1e12;
    end
end

function y = gaussian_1d(p, x)
    sigma = max(abs(p(3)), eps);
    y = p(1) * exp(-((x - p(2)).^2) ./ (2 * sigma^2)) + p(4);
end

function plot_gridmask_profiles(profile, feature_table, cfg, profile_dir, file_tag, algo_name)
    fig = figure('Visible', 'off', 'Position', [100, 100, 1100, 450], 'Color', 'w');

    subplot(1, 2, 1);
    plot(cfg.grid.bin_centers, profile.raw_x, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    hold on;
    plot(cfg.grid.bin_centers, profile.smooth_x, 'b-', 'LineWidth', 1.8);
    plot_feature_markers(feature_table, 'X');
    xlabel('X (mm)');
    ylabel(profile.value_label, 'Interpreter', 'none');
    grid on;
    title(sprintf('X profile at Y = %.2f mm', profile.peak_y), 'Interpreter', 'none');

    subplot(1, 2, 2);
    plot(cfg.grid.bin_centers, profile.raw_y, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1);
    hold on;
    plot(cfg.grid.bin_centers, profile.smooth_y, 'r-', 'LineWidth', 1.8);
    plot_feature_markers(feature_table, 'Y');
    xlabel('Y (mm)');
    ylabel(profile.value_label, 'Interpreter', 'none');
    grid on;
    title(sprintf('Y profile at X = %.2f mm', profile.peak_x), 'Interpreter', 'none');

    exportgraphics(fig, fullfile(profile_dir, ...
        sprintf('GridmaskProfile_%s_%s.png', algo_name, file_tag)), 'Resolution', 300);
    close(fig);
end

function plot_feature_markers(feature_table, axis_name)
    if isempty(feature_table)
        return;
    end
    axis_mask = strcmp(feature_table.Axis, axis_name);
    peak_mask = axis_mask & strcmp(feature_table.FeatureType, 'peak');
    valley_mask = axis_mask & strcmp(feature_table.FeatureType, 'valley');

    if any(peak_mask)
        plot(feature_table.Position_mm(peak_mask), feature_table.Value(peak_mask), ...
            'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 5);
    end
    if any(valley_mask)
        plot(feature_table.Position_mm(valley_mask), feature_table.Value(valley_mask), ...
            'kv', 'MarkerFaceColor', 'c', 'MarkerSize', 5);
    end
end

%% ========================================================================
% Slit post-processing
% =========================================================================
function postprocess_slit(localized, recon, cfg, dirs, file_idx, file_tag, algo_name)
    slit_dir = fullfile(dirs.slit, 'aligned_profiles');
    true_x = localized.true_x_mm(:);
    true_y = localized.true_y_mm(:);
    event_weights = localized.event_weights(:);

    y_edges = -cfg.slit.y_slice_half_range_mm : cfg.slit.y_slice_step_mm : cfg.slit.y_slice_half_range_mm;
    num_slices = numel(y_edges) - 1;
    slice_centers = zeros(num_slices, 1);
    detected_peaks = NaN(num_slices, 1);

    fine_x_edges = -cfg.slit.peak_hist_x_half_range_mm : ...
        cfg.slit.peak_hist_bin_width_mm : cfg.slit.peak_hist_x_half_range_mm;
    if fine_x_edges(end) < cfg.slit.peak_hist_x_half_range_mm
        fine_x_edges = [fine_x_edges, cfg.slit.peak_hist_x_half_range_mm];
    end
    fine_x_centers = (fine_x_edges(1:end-1) + fine_x_edges(2:end)) / 2;

    aligned_x_all = [];
    aligned_w_all = [];

    for ss = 1:num_slices
        slice_centers(ss) = (y_edges(ss) + y_edges(ss+1)) / 2;
        mask = true_y >= y_edges(ss) & true_y < y_edges(ss+1) & abs(true_x) <= cfg.grid.axis_limit;
        local_x = true_x(mask);
        local_w = event_weights(mask);

        if isempty(local_x) || sum(local_w) < cfg.slit.min_slice_weight
            continue;
        end

        counts = weighted_histcounts(local_x, local_w, fine_x_edges);
        counts_smooth = smooth_gaussian_1d(counts, cfg.slit.peak_smooth_sigma_bins);
        [~, max_idx] = max(counts_smooth);
        rough_peak = fine_x_centers(max_idx);

        roi = abs(local_x - rough_peak) < cfg.slit.peak_refine_half_width_mm;
        if any(roi) && sum(local_w(roi)) > 0
            fine_peak = sum(local_x(roi) .* local_w(roi)) / sum(local_w(roi));
        else
            fine_peak = rough_peak;
        end

        detected_peaks(ss) = fine_peak;
        aligned_x_all = [aligned_x_all; local_x - fine_peak]; %#ok<AGROW>
        aligned_w_all = [aligned_w_all; local_w]; %#ok<AGROW>
    end

    eval_x_centers = -cfg.slit.eval_x_half_range_mm : cfg.slit.eval_x_bin_width_mm : cfg.slit.eval_x_half_range_mm;
    eval_x_edges = [eval_x_centers - cfg.slit.eval_x_bin_width_mm / 2, ...
        eval_x_centers(end) + cfg.slit.eval_x_bin_width_mm / 2];
    aligned_counts = weighted_histcounts(aligned_x_all, aligned_w_all, eval_x_edges);

    p0 = [max(aligned_counts), 0, cfg.slit.initial_sigma_mm, min(aligned_counts)];
    p_opt = fit_gaussian_1d_counts(eval_x_centers, aligned_counts, p0);
    fwhm_mm = 2.3548 * abs(p_opt(3));

    fig = figure('Visible', 'off', 'Position', [100, 150, 1200, 500], 'Color', 'w');
    ax1 = subplot(1, 2, 1);
    imagesc(cfg.grid.bin_centers, cfg.grid.bin_centers, recon.display_counts);
    set(gca, 'YDir', 'normal');
    colormap(jet(256));
    axis equal tight;
    hold on;
    cb = colorbar;
    ylabel(cb, cfg.render.colorbar_label, 'Interpreter', 'none');
    apply_image_clim(recon.display_counts, cfg.render);
    valid_peak_mask = ~isnan(detected_peaks);
    plot(detected_peaks(valid_peak_mask), slice_centers(valid_peak_mask), 'r-', 'LineWidth', 2);
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('Detected slit centerline', 'Interpreter', 'none');

    ax2 = subplot(1, 2, 2);
    h_bar = bar(eval_x_centers, aligned_counts, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    hold on;
    h_fit = plot(eval_x_centers, gaussian_1d(p_opt, eval_x_centers), 'r-', 'LineWidth', 2);
    xlabel('Relative distance to center (mm)');
    ylabel('Weighted counts');
    grid on;
    axis square;
    legend([h_bar, h_fit], {'Aligned counts', 'Gaussian fit'}, 'Location', 'northeast');
    xl = xlim;
    yl = ylim;
    text(xl(1) + 0.05 * diff(xl), yl(2) - 0.12 * diff(yl), ...
        sprintf('FWHM = %.2f mm', fwhm_mm), ...
        'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 4);

    ax1.Position = [0.08, 0.16, 0.31, 0.68];
    cb.Position = [0.41, 0.16, 0.014, 0.68];
    ax2.Position = [0.50, 0.16, 0.34, 0.68];

    save_path = fullfile(slit_dir, sprintf('SlitAlignedProfile_%s_%s.png', algo_name, file_tag));
    exportgraphics(fig, save_path, 'Resolution', 300);
    close(fig);

    SlitSummary = table(file_idx, fwhm_mm, p_opt(1), p_opt(2), abs(p_opt(3)), p_opt(4), ...
        'VariableNames', {'FileIndex', 'FWHM_mm', 'FitAmplitude', 'FitCenter_mm', 'FitSigma_mm', 'FitBackground'});
    SliceTable = table(slice_centers, detected_peaks, 'VariableNames', {'SliceCenterY_mm', 'DetectedPeakX_mm'});
    ProfileTable = table(eval_x_centers(:), aligned_counts(:), gaussian_1d(p_opt, eval_x_centers(:)), ...
        'VariableNames', {'RelativeX_mm', 'AlignedCounts', 'GaussianFit'});

    writetable(SlitSummary, fullfile(slit_dir, sprintf('SlitSummary_%s_%s.csv', algo_name, file_tag)));
    writetable(SliceTable, fullfile(slit_dir, sprintf('SlitSliceCenters_%s_%s.csv', algo_name, file_tag)));
    writetable(ProfileTable, fullfile(slit_dir, sprintf('SlitAlignedProfile_%s_%s.csv', algo_name, file_tag)));
    save(fullfile(slit_dir, sprintf('SlitAlignedProfile_%s_%s.mat', algo_name, file_tag)), ...
        'SlitSummary', 'SliceTable', 'ProfileTable', 'aligned_x_all', 'aligned_w_all', 'p_opt');
end

function counts = weighted_histcounts(values, weights, edges)
    counts = zeros(1, numel(edges) - 1);
    if isempty(values)
        return;
    end
    bin_idx = discretize(values, edges);
    valid = ~isnan(bin_idx);
    if any(valid)
        counts = accumarray(bin_idx(valid), weights(valid), [numel(edges)-1, 1])';
    end
end

function smoothed = smooth_gaussian_1d(values, sigma_bins)
    if sigma_bins <= 0
        smoothed = values;
        return;
    end
    radius_bins = max(1, ceil(4 * sigma_bins));
    x = -radius_bins:radius_bins;
    kernel = exp(-(x.^2) ./ (2 * sigma_bins^2));
    kernel = kernel ./ (sum(kernel) + eps);
    smoothed = conv(double(values), kernel, 'same');
end

function p_opt = fit_gaussian_1d_counts(x, y, p0)
    obj = @(p) gaussian_sse_with_penalty(p, x, y, min(x), max(x));
    options = optimset('Display', 'off');
    p_opt = fminsearch(obj, p0, options);
end

%% ========================================================================
% Single-hole post-processing
% =========================================================================
function state = init_single_hole_state(cfg)
    state.composite = zeros(cfg.display_resolution, cfg.display_resolution);
    state.num_images = 0;
    state.file_indices = [];
    state.file_tags = {};
end

function state = accumulate_single_hole(state, recon, cfg, dirs, file_idx, file_tag, algo_name)
    img = recon.display_counts;
    normalized_img = normalize_image_for_composite(img, cfg.single_hole.normalization_method);

    state.composite = state.composite + normalized_img;
    state.num_images = state.num_images + 1;
    state.file_indices(end+1) = file_idx;
    state.file_tags{end+1} = file_tag;

    per_file_path = fullfile(dirs.single_hole, 'per_file', ...
        sprintf('SingleHoleNormalized_%s_%s.png', algo_name, file_tag));
    plot_reconstruction_image(normalized_img, cfg.grid.bin_centers, cfg.grid.axis_limit, ...
        per_file_path, single_hole_render_cfg(cfg.render));
end

function img = normalize_image_for_composite(img, method)
    switch lower(method)
        case 'peak'
            denom = max(img(:));
        case 'integral'
            denom = sum(img(:));
        otherwise
            error('Unknown single-hole normalization method: %s', method);
    end

    if denom > 0
        img = img ./ denom;
    end
end

function render_cfg = single_hole_render_cfg(render_cfg)
    render_cfg.colorbar_label = 'Normalized intensity';
    render_cfg.colorbar_mode = 'normalized';
end

function finalize_single_hole(state, cfg, dirs, algo_name)
    if isempty(state) || state.num_images == 0
        warning('No single-hole images were accumulated.');
        return;
    end

    composite_path = fullfile(dirs.single_hole, 'composite', ...
        sprintf('SingleHoleComposite_%s_N%d.png', algo_name, state.num_images));
    plot_reconstruction_image(state.composite, cfg.grid.bin_centers, cfg.grid.axis_limit, ...
        composite_path, single_hole_render_cfg(cfg.render));

    save(fullfile(dirs.single_hole, 'composite', ...
        sprintf('SingleHoleComposite_%s_N%d.mat', algo_name, state.num_images)), ...
        'state', 'cfg');
end

%% ========================================================================
% Single-hole scan position-linearity analysis
% =========================================================================
function state = init_single_hole_scan_state(cfg)
    if numel(cfg.single_hole_scan.file_indices) ~= numel(cfg.single_hole_scan.mechanical_y_mm)
        error('single_hole_scan.file_indices and single_hole_scan.mechanical_y_mm must have the same length.');
    end
    state.composite = zeros(cfg.display_resolution, cfg.display_resolution);
    state.num_images = 0;
    state.rows = cell(0, 16);
    state.file_tags = {};
end

function state = accumulate_single_hole_scan(state, localized, recon, cfg, dirs, file_idx, file_tag, algo_name)
    scan_pos = find(cfg.single_hole_scan.file_indices == file_idx, 1);
    if isempty(scan_pos)
        warning('No mechanical position is configured for single-hole file p%04d.', file_idx);
        return;
    end

    mechanical_x = cfg.single_hole_scan.mechanical_x_mm;
    mechanical_y = cfg.single_hole_scan.mechanical_y_mm(scan_pos);
    source_name = lower(strtrim(char(cfg.single_hole_scan.fit_image_source)));
    if ~isfield(recon, source_name)
        error('Unknown single-hole scan fit image source: %s', source_name);
    end
    fit = fit_single_hole_spot_2d(recon.(source_name), cfg);

    img = recon.display_counts;
    normalized_img = normalize_image_for_composite(img, cfg.single_hole.normalization_method);
    state.composite = state.composite + normalized_img;
    state.num_images = state.num_images + 1;
    state.file_tags{end+1} = file_tag;

    error_x = fit.center_x_mm - mechanical_x;
    error_y = fit.center_y_mm - mechanical_y;
    state.rows(end+1, :) = {file_idx, mechanical_x, mechanical_y, fit.center_x_mm, ...
        fit.center_y_mm, error_x, error_y, abs(error_y), hypot(error_x, error_y), ...
        fit.sigma_x_mm, fit.sigma_y_mm, fit.amplitude, fit.background, fit.r_squared, ...
        fit.valid, numel(localized.true_x_mm)}; %#ok<AGROW>

    per_file_path = fullfile(dirs.single_hole_scan, 'per_file', ...
        sprintf('SingleHoleScan_%s_%s.png', algo_name, file_tag));
    plot_single_hole_scan_frame(normalized_img, fit, mechanical_x, mechanical_y, ...
        cfg, per_file_path);
end

function fit = fit_single_hole_spot_2d(img_data, cfg)
    img = double(img_data);
    img(~isfinite(img)) = 0;
    if isempty(img) || max(img(:)) <= 0
        fit = empty_single_hole_fit();
        return;
    end

    [~, max_idx] = max(img(:));
    [seed_y_idx, seed_x_idx] = ind2sub(size(img), max_idx);
    seed_x = cfg.grid.bin_centers(seed_x_idx);
    seed_y = cfg.grid.bin_centers(seed_y_idx);
    half_width = cfg.single_hole_scan.fit_roi_half_width_mm;
    x_mask = abs(cfg.grid.bin_centers - seed_x) <= half_width;
    y_mask = abs(cfg.grid.bin_centers - seed_y) <= half_width;
    [X, Y] = meshgrid(cfg.grid.bin_centers(x_mask), cfg.grid.bin_centers(y_mask));
    Z = img(y_mask, x_mask);
    z_scale = max(Z(:));
    if isempty(Z) || z_scale <= 0
        fit = empty_single_hole_fit();
        return;
    end
    Z = Z ./ z_scale;

    background0 = min(Z(:));
    amplitude0 = max(Z(:)) - background0;
    p0 = [amplitude0, seed_x, seed_y, cfg.single_hole_scan.initial_sigma_mm, ...
        cfg.single_hole_scan.initial_sigma_mm, background0];
    bounds = [min(X(:)), max(X(:)), min(Y(:)), max(Y(:))];
    objective = @(p) single_hole_gaussian_sse(p, X, Y, Z, cfg.single_hole_scan, bounds);
    options = optimset('Display', 'off', 'MaxIter', 500, 'MaxFunEvals', 2000);
    p_opt = fminsearch(objective, p0, options);

    model = single_hole_gaussian_2d(p_opt, X, Y);
    residual_sse = sum((Z(:) - model(:)).^2);
    total_sse = sum((Z(:) - mean(Z(:))).^2);
    r_squared = 1 - residual_sse / (total_sse + eps);
    sigma_x = abs(p_opt(4));
    sigma_y = abs(p_opt(5));
    valid = p_opt(1) >= 0 && p_opt(6) >= 0 && ...
        p_opt(2) >= bounds(1) && p_opt(2) <= bounds(2) && ...
        p_opt(3) >= bounds(3) && p_opt(3) <= bounds(4) && ...
        sigma_x >= cfg.single_hole_scan.min_sigma_mm && ...
        sigma_x <= cfg.single_hole_scan.max_sigma_mm && ...
        sigma_y >= cfg.single_hole_scan.min_sigma_mm && ...
        sigma_y <= cfg.single_hole_scan.max_sigma_mm;

    fit.center_x_mm = p_opt(2);
    fit.center_y_mm = p_opt(3);
    fit.sigma_x_mm = sigma_x;
    fit.sigma_y_mm = sigma_y;
    fit.amplitude = p_opt(1) * z_scale;
    fit.background = p_opt(6) * z_scale;
    fit.r_squared = r_squared;
    fit.valid = valid;
end

function fit = empty_single_hole_fit()
    fit.center_x_mm = NaN;
    fit.center_y_mm = NaN;
    fit.sigma_x_mm = NaN;
    fit.sigma_y_mm = NaN;
    fit.amplitude = NaN;
    fit.background = NaN;
    fit.r_squared = NaN;
    fit.valid = false;
end

function err = single_hole_gaussian_sse(p, X, Y, Z, scan_cfg, bounds)
    model = single_hole_gaussian_2d(p, X, Y);
    err = sum((Z(:) - model(:)).^2);
    sigma_x = abs(p(4));
    sigma_y = abs(p(5));
    invalid = p(1) < 0 || p(6) < 0 || ...
        p(2) < bounds(1) || p(2) > bounds(2) || ...
        p(3) < bounds(3) || p(3) > bounds(4) || ...
        sigma_x < scan_cfg.min_sigma_mm || sigma_x > scan_cfg.max_sigma_mm || ...
        sigma_y < scan_cfg.min_sigma_mm || sigma_y > scan_cfg.max_sigma_mm;
    if invalid
        err = err + 1e6;
    end
end

function model = single_hole_gaussian_2d(p, X, Y)
    sigma_x = max(abs(p(4)), eps);
    sigma_y = max(abs(p(5)), eps);
    model = p(6) + p(1) .* exp(-((X - p(2)).^2 ./ (2 * sigma_x^2) + ...
        (Y - p(3)).^2 ./ (2 * sigma_y^2)));
end

function plot_single_hole_scan_frame(img, fit, mechanical_x, mechanical_y, cfg, save_path)
    fig = figure('Visible', 'off', 'Position', [100, 100, 650, 600], 'Color', 'w');
    imagesc(cfg.grid.bin_centers, cfg.grid.bin_centers, img);
    set(gca, 'YDir', 'normal', 'XTick', -25:5:25, 'YTick', -25:5:25, 'Layer', 'top');
    axis equal tight;
    colormap(jet(256));
    cb = colorbar;
    ylabel(cb, 'Normalized intensity', 'Interpreter', 'none');
    hold on;
    plot(mechanical_x, mechanical_y, 'wo', 'LineWidth', 1.3, 'MarkerSize', 9);
    if fit.valid
        plot(fit.center_x_mm, fit.center_y_mm, 'rx', 'LineWidth', 1.8, 'MarkerSize', 10);
    end
    xlabel('X (mm)');
    ylabel('Y (mm)');
    exportgraphics(fig, save_path, 'Resolution', 300);
    close(fig);
end

function finalize_single_hole_scan(state, cfg, dirs, algo_name)
    if isempty(state) || state.num_images == 0 || isempty(state.rows)
        warning('No single-hole scan images were accumulated.');
        return;
    end

    names = {'FileIndex', 'MechanicalX_mm', 'MechanicalY_mm', 'FittedX_mm', ...
        'FittedY_mm', 'ErrorX_mm', 'ErrorY_mm', 'AbsoluteYError_mm', ...
        'EuclideanError_mm', 'SigmaX_mm', 'SigmaY_mm', 'Amplitude', ...
        'Background', 'FitR2', 'FitValid', 'NumValidEvents'};
    PositionTable = cell2table(state.rows, 'VariableNames', names);
    PositionTable.FitValid = logical(PositionTable.FitValid) & ...
        isfinite(PositionTable.FittedY_mm) & isfinite(PositionTable.FittedX_mm);
    fwhm_factor = 2 * sqrt(2 * log(2));
    PositionTable.FWHM_X_mm = fwhm_factor .* PositionTable.SigmaX_mm;
    PositionTable.FWHM_Y_mm = fwhm_factor .* PositionTable.SigmaY_mm;
    PositionTable.FWHM_Mean_mm = (PositionTable.FWHM_X_mm + PositionTable.FWHM_Y_mm) ./ 2;
    PositionTable.FWHM_GeometricMean_mm = sqrt(PositionTable.FWHM_X_mm .* PositionTable.FWHM_Y_mm);
    PositionTable = sortrows(PositionTable, 'MechanicalY_mm');
    valid = PositionTable.FitValid;
    if nnz(valid) < 2
        warning('Fewer than two valid single-hole fits were available for %s.', algo_name);
        return;
    end

    mechanical_y = PositionTable.MechanicalY_mm(valid);
    fitted_y = PositionTable.FittedY_mm(valid);
    metrics = compute_single_hole_linearity_metrics(mechanical_y, fitted_y);
    coefficients = [metrics.slope, metrics.intercept_mm];
    PositionTable.FittedLineY_mm = polyval(coefficients, PositionTable.MechanicalY_mm);
    PositionTable.LineFitResidualY_mm = PositionTable.FittedY_mm - PositionTable.FittedLineY_mm;
    PositionTable.OffsetCorrectedYError_mm = NaN(height(PositionTable), 1);
    PositionTable.OffsetCorrectedYError_mm(valid) = ...
        PositionTable.ErrorY_mm(valid) - metrics.mean_y_error_mm;
    r_squared = metrics.r_squared;
    mae_mm = mean(abs(PositionTable.ErrorY_mm(valid)));
    rmse_mm = sqrt(mean(PositionTable.ErrorY_mm(valid).^2));
    max_abs_error_mm = max(abs(PositionTable.ErrorY_mm(valid)));
    mean_x = mean(PositionTable.FittedX_mm(valid));

    LinearitySummary = table({algo_name}, nnz(valid), coefficients(1), coefficients(2), ...
        metrics.slope_std_error, metrics.intercept_std_error_mm, ...
        r_squared, mae_mm, rmse_mm, max_abs_error_mm, ...
        metrics.mean_y_error_mm, metrics.offset_corrected_error_sd_mm, ...
        metrics.line_fit_residual_sd_mm, mean_x, ...
        mean_x - cfg.single_hole_scan.mechanical_x_mm, ...
        'VariableNames', {'Algorithm', 'NumPositions', 'Slope', 'Intercept_mm', ...
        'SlopeStdError', 'InterceptStdError_mm', 'R2', 'MAE_Y_mm', ...
        'RMSE_Y_mm', 'MaxAbsoluteError_Y_mm', 'MeanYOffset_mm', ...
        'OffsetCorrectedErrorSD_mm', 'LineFitResidualSD_mm', ...
        'MeanFittedX_mm', 'MeanXOffset_mm'});
    Table2Summary = make_single_hole_table2_summary({algo_name}, LinearitySummary);

    valid_sr = PositionTable.FitValid & isfinite(PositionTable.FWHM_Mean_mm);
    SpatialResolutionTable = PositionTable(:, {'FileIndex', 'MechanicalX_mm', ...
        'MechanicalY_mm', 'FittedX_mm', 'FittedY_mm', 'FWHM_X_mm', ...
        'FWHM_Y_mm', 'FWHM_Mean_mm', 'FWHM_GeometricMean_mm', ...
        'FitR2', 'FitValid', 'NumValidEvents'});
    SpatialResolutionTable.Algorithm = repmat({algo_name}, height(SpatialResolutionTable), 1);
    SpatialResolutionSummary = table({algo_name}, nnz(valid_sr), ...
        mean(PositionTable.FWHM_X_mm(valid_sr)), ...
        mean(PositionTable.FWHM_Y_mm(valid_sr)), ...
        mean(PositionTable.FWHM_Mean_mm(valid_sr)), ...
        std(PositionTable.FWHM_Mean_mm(valid_sr)), ...
        min(PositionTable.FWHM_Mean_mm(valid_sr)), ...
        max(PositionTable.FWHM_Mean_mm(valid_sr)), ...
        'VariableNames', {'Algorithm', 'NumPositions', 'MeanFWHM_X_mm', ...
        'MeanFWHM_Y_mm', 'MeanFWHM_Mean_mm', 'StdFWHM_Mean_mm', ...
        'MinFWHM_Mean_mm', 'MaxFWHM_Mean_mm'});

    table_dir = fullfile(dirs.single_hole_scan, 'tables');
    figure_dir = fullfile(dirs.single_hole_scan, 'figures');
    if isfield(state, 'file_tags') && ~isempty(state.file_tags)
        tag = sprintf('%s_%s', algo_name, summarize_file_tags(state.file_tags));
    else
        tag = sprintf('%s_p%04d_p%04d', algo_name, min(PositionTable.FileIndex), max(PositionTable.FileIndex));
    end
    writetable(PositionTable, fullfile(table_dir, sprintf('SingleHoleScanPositions_%s.csv', tag)));
    writetable(LinearitySummary, fullfile(table_dir, sprintf('SingleHoleScanSummary_%s.csv', tag)));
    writetable(Table2Summary, fullfile(table_dir, sprintf('SingleHoleScanTable2_%s.csv', tag)));
    writetable(SpatialResolutionTable, fullfile(table_dir, ...
        sprintf('SingleHoleScanSpatialResolution_%s.csv', tag)));
    writetable(SpatialResolutionSummary, fullfile(table_dir, ...
        sprintf('SingleHoleScanSpatialResolutionSummary_%s.csv', tag)));

    composite_img = state.composite ./ (max(state.composite(:)) + eps);
    plot_single_hole_scan_linearity_figure(composite_img, PositionTable, LinearitySummary, ...
        cfg, fullfile(figure_dir, sprintf('SingleHoleScanLinearity_%s.png', tag)));
    save(fullfile(table_dir, sprintf('SingleHoleScanData_%s.mat', tag)), ...
        'PositionTable', 'LinearitySummary', 'SpatialResolutionTable', ...
        'SpatialResolutionSummary', 'Table2Summary', 'state', 'cfg');

    plot_single_hole_scan_comparison_if_available(cfg, dirs);
end

function metrics = compute_single_hole_linearity_metrics(mechanical_y, fitted_y)
    x = double(mechanical_y(:));
    y = double(fitted_y(:));
    n = numel(x);

    X = [x, ones(n, 1)];
    beta = X \ y;
    model_y = X * beta;
    residual_to_line = y - model_y;
    dof = n - 2;

    if dof > 0
        mse = sum(residual_to_line.^2) / dof;
        cov_beta = mse * ((X' * X) \ eye(2));
        slope_std_error = sqrt(cov_beta(1, 1));
        intercept_std_error_mm = sqrt(cov_beta(2, 2));
        line_fit_residual_sd_mm = sqrt(mse);
    else
        slope_std_error = NaN;
        intercept_std_error_mm = NaN;
        line_fit_residual_sd_mm = NaN;
    end

    y_error = y - x;
    mean_y_error_mm = mean(y_error);
    if n > 1
        offset_corrected_error_sd_mm = std(y_error - mean_y_error_mm, 0);
    else
        offset_corrected_error_sd_mm = NaN;
    end

    metrics.slope = beta(1);
    metrics.intercept_mm = beta(2);
    metrics.slope_std_error = slope_std_error;
    metrics.intercept_std_error_mm = intercept_std_error_mm;
    metrics.r_squared = 1 - sum(residual_to_line.^2) / (sum((y - mean(y)).^2) + eps);
    metrics.mean_y_error_mm = mean_y_error_mm;
    metrics.offset_corrected_error_sd_mm = offset_corrected_error_sd_mm;
    metrics.line_fit_residual_sd_mm = line_fit_residual_sd_mm;
end

function Table2Summary = make_single_hole_table2_summary(method_names, linearity_summary)
    Table2Summary = table(method_names(:), linearity_summary.Slope, ...
        linearity_summary.SlopeStdError, linearity_summary.Intercept_mm, ...
        linearity_summary.InterceptStdError_mm, ...
        linearity_summary.OffsetCorrectedErrorSD_mm, ...
        'VariableNames', {'Method', 'k', 'Delta_k', 'b_mm', 'Delta_b_mm', 'SD_mm'});
end

function plot_single_hole_scan_linearity_figure(composite_img, positions, summary, cfg, save_path)
    valid = positions.FitValid;
    fig = figure('Visible', 'off', 'Position', [100, 100, 1180, 500], 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    imagesc(cfg.grid.bin_centers, cfg.grid.bin_centers, composite_img);
    set(gca, 'YDir', 'normal', 'XTick', -25:5:25, 'YTick', -25:5:25, 'Layer', 'top');
    axis equal tight;
    colormap(jet(256));
    hold on;
    plot(repmat(cfg.single_hole_scan.mechanical_x_mm, numel(cfg.single_hole_scan.mechanical_y_mm), 1), ...
        cfg.single_hole_scan.mechanical_y_mm(:), 'w--', 'LineWidth', 1.3);
    plot(positions.FittedX_mm(valid), positions.FittedY_mm(valid), 'ro-', ...
        'LineWidth', 1.3, 'MarkerFaceColor', 'r', 'MarkerSize', 5);
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('Reconstructed scan path');

    nexttile;
    hold on;
    ideal_y = cfg.single_hole_scan.mechanical_y_mm(:);
    plot(ideal_y, ideal_y, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Ideal');
    plot(positions.MechanicalY_mm(valid), positions.FittedY_mm(valid), 'o', ...
        'Color', [0.00, 0.45, 0.74], 'MarkerFaceColor', [0.00, 0.45, 0.74], ...
        'MarkerSize', 6, 'DisplayName', 'Gaussian center');
    plot(positions.MechanicalY_mm(valid), positions.FittedLineY_mm(valid), '-', ...
        'Color', [0.85, 0.33, 0.10], 'LineWidth', 1.6, 'DisplayName', 'Linear fit');
    xlabel('Mechanical Y position (mm)');
    ylabel('Reconstructed Y position (mm)');
    xlim([-22, 22]);
    ylim([-25, 25]);
    axis square;
    grid on;
    legend('Location', 'northwest');
    text(-20, -16, sprintf('slope = %.2f\nb = %.2f mm\nR^2 = %.4f\nMAE = %.3f mm', ...
        summary.Slope, summary.Intercept_mm, summary.R2, summary.MAE_Y_mm), ...
        'BackgroundColor', 'w', 'Margin', 4, 'FontSize', 10);
    title('Position linearity');

    exportgraphics(fig, save_path, 'Resolution', 300);
    close(fig);
end

function plot_single_hole_scan_comparison_if_available(cfg, dirs)
    tags = {'Anger_Standard', 'Fabbri_AnalyticalLSF', 'LSE_Gaussian', 'LSE_Gaussian_64ch'};
    labels = {'Anger', 'Analytical least-square fitting', 'Proposed method', ...
        '64-channel LSE softmax'};
    table2_labels = {'Anger centroid', 'Analytical least-square fitting', ...
        'Proposed MC-library-based method', '64-channel LSE softmax'};
    colors = {[0.25, 0.25, 0.25], [0.85, 0.33, 0.10], [0.00, 0.45, 0.74], ...
        [0.47, 0.67, 0.19]};
    table_dir = fullfile(dirs.single_hole_scan, 'tables');
    all_positions = table();
    all_summary = table();
    all_table2 = table();
    available = false(1, numel(tags));

    for tt = 1:numel(tags)
        position_path = fullfile(table_dir, sprintf('SingleHoleScanPositions_%s_p0801_p0809.csv', tags{tt}));
        summary_path = fullfile(table_dir, sprintf('SingleHoleScanSummary_%s_p0801_p0809.csv', tags{tt}));
        if exist(position_path, 'file') && exist(summary_path, 'file')
            positions = readtable(position_path);
            summary = readtable(summary_path);
            positions = positions(:, {'FileIndex', 'MechanicalY_mm', ...
                'FittedY_mm', 'ErrorY_mm', 'FitValid'});
            positions.Algorithm = repmat(labels(tt), height(positions), 1);
            summary.DisplayName = repmat(labels(tt), height(summary), 1);
            all_positions = [all_positions; positions]; %#ok<AGROW>
            all_summary = append_table_with_union(all_summary, summary);
            valid = logical(positions.FitValid) & ...
                isfinite(positions.MechanicalY_mm) & isfinite(positions.FittedY_mm);
            if nnz(valid) >= 2
                metrics = compute_single_hole_linearity_metrics( ...
                    positions.MechanicalY_mm(valid), positions.FittedY_mm(valid));
                method_table2 = table(table2_labels(tt), metrics.slope, ...
                    metrics.slope_std_error, metrics.intercept_mm, ...
                    metrics.intercept_std_error_mm, metrics.offset_corrected_error_sd_mm, ...
                    'VariableNames', {'Method', 'k', 'Delta_k', 'b_mm', 'Delta_b_mm', 'SD_mm'});
                all_table2 = [all_table2; method_table2]; %#ok<AGROW>
            end
            available(tt) = true;
        end
    end
    if ~any(available)
        return;
    end

    fig = figure('Visible', 'off', 'Position', [100, 100, 1000, 450], 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_linearity = nexttile;
    hold(ax_linearity, 'on');
    ideal_y = cfg.single_hole_scan.mechanical_y_mm(:);
    plot(ax_linearity, ideal_y, ideal_y, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Ideal');
    ax_error = nexttile;
    hold(ax_error, 'on');

    for tt = 1:numel(tags)
        if ~available(tt)
            continue;
        end
        mask = strcmp(all_positions.Algorithm, labels{tt}) & logical(all_positions.FitValid);
        positions = all_positions(mask, :);
        plot(ax_linearity, positions.MechanicalY_mm, positions.FittedY_mm, 'o-', ...
            'Color', colors{tt}, 'MarkerFaceColor', colors{tt}, ...
            'LineWidth', 1.4, 'DisplayName', labels{tt});
        plot(ax_error, positions.MechanicalY_mm, abs(positions.ErrorY_mm), 'o-', ...
            'Color', colors{tt}, 'MarkerFaceColor', colors{tt}, ...
            'LineWidth', 1.4, 'DisplayName', labels{tt});
    end

    xlabel(ax_linearity, 'Mechanical Y position (mm)');
    ylabel(ax_linearity, 'Reconstructed Y position (mm)');
    title(ax_linearity, 'Position linearity');
    grid(ax_linearity, 'on');
    axis(ax_linearity, 'square');
    xlim(ax_linearity, [-22, 22]);
    ylim(ax_linearity, [-25, 25]);
    legend(ax_linearity, 'Location', 'northwest');
    xlabel(ax_error, 'Mechanical Y position (mm)');
    ylabel(ax_error, '|Position error| (mm)');
    title(ax_error, 'Absolute position error');
    grid(ax_error, 'on');
    axis(ax_error, 'square');
    xlim(ax_error, [-22, 22]);
    legend(ax_error, 'Location', 'best');

    figure_path = fullfile(dirs.single_hole_scan, 'figures', ...
        'SingleHoleScanComparison_AvailableMethods.png');
    exportgraphics(fig, figure_path, 'Resolution', 300);
    close(fig);
    writetable(all_summary, fullfile(table_dir, 'SingleHoleScanComparisonSummary_AvailableMethods.csv'));
    if ~isempty(all_table2)
        writetable(all_table2, fullfile(table_dir, 'SingleHoleScanTable2Summary_AvailableMethods.csv'));
    end
end

function out = append_table_with_union(out, next)
    if isempty(out)
        out = next;
        return;
    end

    out_names = out.Properties.VariableNames;
    next_names = next.Properties.VariableNames;
    all_names = unique([out_names, next_names], 'stable');

    out = add_missing_table_columns(out, all_names, next);
    next = add_missing_table_columns(next, all_names, out);
    out = [out(:, all_names); next(:, all_names)];
end

function tbl = add_missing_table_columns(tbl, all_names, template_tbl)
    for ii = 1:numel(all_names)
        name = all_names{ii};
        if ismember(name, tbl.Properties.VariableNames)
            continue;
        end
        template = template_tbl.(name);
        tbl.(name) = missing_values_like(template, height(tbl));
    end
end

function values = missing_values_like(template, n)
    if iscell(template)
        values = repmat({''}, n, 1);
    elseif isstring(template)
        values = strings(n, 1);
    elseif islogical(template)
        values = false(n, 1);
    elseif isnumeric(template)
        values = NaN(n, size(template, 2));
    elseif iscategorical(template)
        values = categorical(repmat({''}, n, 1));
    else
        values = repmat(missing, n, 1);
    end
end
