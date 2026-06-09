function locator = prepare_fabbri_analytical_lsf_locator(cfg)
%PREPARE_FABBRI_ANALYTICAL_LSF_LOCATOR Configure analytical Fabbri fitting.

method_dir = fileparts(mfilename('fullpath'));
fabbri_dir = fullfile(method_dir, 'fabbri_analytical_lsf');
if ~exist(fabbri_dir, 'dir')
    legacy_fabbri_dir = fullfile(cfg.repo_root, 'Fabbri');
    if exist(legacy_fabbri_dir, 'dir')
        fabbri_dir = legacy_fabbri_dir;
    else
        error('Fabbri analytical fitting functions not found: %s', fabbri_dir);
    end
end
addpath(fabbri_dir, '-begin');

fabbri_cfg = fabbri_default_config(cfg.repo_root);
fabbri_cfg.input_data_folder = cfg.input_data_folder;
fabbri_cfg.display_resolution = cfg.display_resolution;
fabbri_cfg.axis_limit_mm = cfg.show_size_mm / 2;
fabbri_cfg.position_bounds_mm = [-cfg.show_size_mm / 2, cfg.show_size_mm / 2];
fabbri_cfg.anode_pitch_mm = cfg.detector_pitch_mm;

if isfield(cfg, 'fabbri')
    fabbri_cfg = merge_struct_fields(fabbri_cfg, cfg.fabbri);
end

locator.method = 'fabbri_analytical_lsf';
locator.algorithm_name = 'Fabbri_AnalyticalLSF';
locator.fabbri_cfg = fabbri_cfg;
end

function out = merge_struct_fields(out, updates)
fields = fieldnames(updates);
for ii = 1:numel(fields)
    out.(fields{ii}) = updates.(fields{ii});
end
end
