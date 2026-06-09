function cfg = IMAGE_default_config(repo_root)
%IMAGE_DEFAULT_CONFIG Build the default cfg used by the IMAGE UI.

if nargin < 1 || isempty(repo_root)
    repo_root = fileparts(mfilename('fullpath'));
end

cfg.script_dir = repo_root;
cfg.repo_root = repo_root;
cfg.input_data_folder = fullfile(repo_root, 'input_data');
cfg.output_root = fullfile(cfg.script_dir, 'IMAGE_Recon_Results');

cfg.calib_root = fullfile(repo_root, 'CALIB');
cfg.lse_calibration_folder = fullfile(cfg.calib_root, 'Calibration_Results_FloodCDF_LSE');
cfg.anger_calibration_folder = fullfile(cfg.calib_root, 'Calibration_Results_AngerCDF');
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

cfg.postprocess_mode = 'gridmask';
cfg.file_basename = 'calib_events_pgridmask%04d.mat';
cfg.file_indices_to_process = 1010;
cfg.selected_files = {};

cfg.enable_cdf_correction = false;
cfg.enable_ucm_correction = false;

cfg.render.enable_gaussian_diffusion = true;
cfg.render.gaussian_sigma_bins = 0.5;
cfg.render.gaussian_radius_bins = 10;
cfg.render.gaussian_kernel_normalization = 'sum1';
cfg.render.colorbar_mode = 'counts_per_mm2';
cfg.render.enable_grid_lines = false;
cfg.render.enable_slit_style_clim = true;
cfg.render.clim_low_percentile = 0;
cfg.render.clim_high_percentile = 99.99;

cfg.localization.method = 'lse_softmax_64ch';
cfg.localization.coarse_radius = 25;
cfg.localization.enable_softmax_interpolation = true;
cfg.localization.lse_temperature = 0.005;
cfg.localization.baseline_ratio = 0.0;
cfg.localization.top_k_ratio = 0.5;
cfg.localization.batch_size = 10000;
cfg.localization.batch_size_64ch = 1000;
cfg.localization.enable_quality_filter = true;
cfg.localization.residual_percentile_cutoff = 80;
cfg.localization.reset_gpu = true;

cfg.anger.position_weight_mode = 'rtp';
cfg.anger.projection_power = 2;
cfg.anger.sensor_centers_mm = [-21, -15, -9, -3, 3, 9, 15, 21];
cfg.anger.clip_negative = true;

cfg.fabbri.max_events_per_file = Inf;
cfg.fabbri.event_stride = 1;
cfg.fabbri.progress_interval = 2000;
cfg.fabbri.shape_b = 1.5;
cfg.fabbri.initial_alpha = 0.015;
cfg.fabbri.max_iterations = 80;
cfg.fabbri.require_convergence = true;
cfg.fabbri.reject_boundary_hits = true;
cfg.fabbri.boundary_tolerance_mm = 1e-6;

cfg.gridmask.manual_peak_search_radius_mm = 3.0;
cfg.gridmask.profile_smooth_window_bins = 20;
cfg.gridmask.min_peak_relative_height = 0.15;
cfg.gridmask.min_peak_distance_mm = 2.0;
cfg.gridmask.fit_half_width_mm = 0.8;
cfg.gridmask.initial_sigma_mm = 0.4;
cfg.gridmask.max_num_peaks = 32;
cfg.gridmask.nominal_spacing_mm = 5.0;
cfg.gridmask.max_adjacent_spacing_factor = 1.5;

cfg.slit.y_slice_half_range_mm = 15.0;
cfg.slit.y_slice_step_mm = 0.1;
cfg.slit.min_slice_weight = 20;
cfg.slit.peak_hist_x_half_range_mm = 3.0;
cfg.slit.peak_hist_bin_width_mm = 0.05;
cfg.slit.peak_smooth_sigma_bins = 1.0;
cfg.slit.peak_refine_half_width_mm = 1.5;
cfg.slit.eval_x_half_range_mm = 15.0;
cfg.slit.eval_x_bin_width_mm = 0.05;
cfg.slit.initial_sigma_mm = 0.8;

cfg.single_hole.normalization_method = 'peak';
cfg.single_hole.composite_method = 'sum_normalized';

cfg.single_hole_scan.file_indices = 801:809;
cfg.single_hole_scan.mechanical_x_mm = 15.0;
cfg.single_hole_scan.mechanical_y_mm = -20:5:20;
cfg.single_hole_scan.fit_image_source = 'raw_counts';
cfg.single_hole_scan.fit_roi_half_width_mm = 3.0;
cfg.single_hole_scan.initial_sigma_mm = 0.8;
cfg.single_hole_scan.min_sigma_mm = 0.05;
cfg.single_hole_scan.max_sigma_mm = 5.0;
end
