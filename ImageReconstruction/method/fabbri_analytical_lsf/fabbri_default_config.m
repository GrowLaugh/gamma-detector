function cfg = fabbri_default_config(repo_root)
%FABBRI_DEFAULT_CONFIG Default parameters for analytical MSE positioning.

if nargin < 1 || isempty(repo_root)
    this_dir = fileparts(mfilename('fullpath'));
    repo_root = fileparts(fileparts(fileparts(this_dir)));
end

cfg.repo_root = repo_root;
cfg.input_data_folder = fullfile(repo_root, 'step1rawDataProcess', 'calib_data1_autofit');
cfg.file_basename = 'calib_events_p%04d.mat';
cfg.file_indices = 503;
cfg.output_folder = fullfile(repo_root, 'IMAGE', 'IMAGE_Recon_Results', ...
    'fabbri_analytical_lsf');

cfg.max_events_per_file = 20000;
cfg.event_stride = 1;
cfg.progress_interval = 2000;

cfg.clip_negative = true;
cfg.normalize_each_event = true;

% H12700/H8500 style 8 x 8 anode geometry. The existing LSE code uses 6 mm,
% so that is the default here for a direct comparison.
cfg.anode_pitch_mm = 6.0;
cfg.position_bounds_mm = [-25, 25];
cfg.axis_limit_mm = 25;
cfg.display_resolution = 512;

% Modified Scrimger-Baker model.
cfg.shape_b = 1.5;
cfg.initial_alpha = 0.015;
cfg.alpha_bounds = [1e-4, 0.3];

% Iterative minimizer.
cfg.max_iterations = 80;
cfg.position_tolerance_mm = 0.01;
cfg.objective_tolerance = 1e-10;
cfg.initial_lambda = 1e-3;
cfg.max_inner_trials = 8;
cfg.max_position_step_mm = 2.0;
cfg.max_log_step = 0.7;

% Event acceptance.
cfg.enable_residual_filter = true;
cfg.residual_percentile_cutoff = 80;
cfg.reject_boundary_hits = false;
cfg.boundary_tolerance_mm = 1e-6;

% Image and figure generation.
cfg.enable_normalization = true;
cfg.enable_grid_lines = false;
cfg.enable_fwhm_analysis = true;
cfg.show_figures = false;
cfg.enable_slit_style_clim = true;
cfg.clim_low_percentile = 0;
cfg.clim_high_percentile = 99.99;
end
