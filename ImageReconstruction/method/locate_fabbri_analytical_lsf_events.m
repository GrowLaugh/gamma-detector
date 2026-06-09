function localized = locate_fabbri_analytical_lsf_events(planeset, locator, ~)
%LOCATE_FABBRI_ANALYTICAL_LSF_EVENTS Fit events with the analytical Fabbri model.

fabbri_cfg = locator.fabbri_cfg;
num_events = numel(planeset);

event_indices = 1:num_events;
if isfield(fabbri_cfg, 'event_stride') && fabbri_cfg.event_stride > 1
    event_indices = event_indices(1:fabbri_cfg.event_stride:end);
end
if isfield(fabbri_cfg, 'max_events_per_file') && isfinite(fabbri_cfg.max_events_per_file)
    event_indices = event_indices(1:min(numel(event_indices), fabbri_cfg.max_events_per_file));
end

planeset_subset = planeset(event_indices);
fprintf('Fabbri analytical fitting: %d/%d events selected.\n', numel(planeset_subset), num_events);
[coords, fitinfo] = fabbri_position_events(planeset_subset, fabbri_cfg);

finite_mask = all(isfinite(coords), 2) & isfinite(fitinfo.residual_mse);
if isfield(fabbri_cfg, 'require_convergence') && fabbri_cfg.require_convergence
    accepted_mask = finite_mask & fitinfo.converged;
else
    accepted_mask = finite_mask;
end
accepted_before_boundary_filter = accepted_mask;

boundary_hit_mask = false(size(accepted_mask));
if isfield(fabbri_cfg, 'reject_boundary_hits') && fabbri_cfg.reject_boundary_hits
    tol = 1e-6;
    if isfield(fabbri_cfg, 'boundary_tolerance_mm') && ~isempty(fabbri_cfg.boundary_tolerance_mm)
        tol = fabbri_cfg.boundary_tolerance_mm;
    end
    bounds = double(fabbri_cfg.position_bounds_mm(:)');
    x = double(coords(:, 1));
    y = double(coords(:, 2));
    boundary_hit_mask = x <= bounds(1) + tol | x >= bounds(2) - tol | ...
        y <= bounds(1) + tol | y >= bounds(2) - tol;
    accepted_mask = accepted_mask & ~boundary_hit_mask;
end

localized.raw_x_mm = double(coords(accepted_mask, 1));
localized.raw_y_mm = double(coords(accepted_mask, 2));
localized.residual = double(fitinfo.residual_mse(accepted_mask));
localized.valid_mask = true(nnz(accepted_mask), 1);
localized.method_info.event_indices = event_indices(:);
localized.method_info.accepted_event_indices = event_indices(accepted_mask);
localized.method_info.fitinfo = fitinfo;
localized.method_info.accepted_mask = accepted_mask;
localized.method_info.boundary_hit_mask = boundary_hit_mask;
localized.method_info.num_boundary_hits_rejected = nnz(accepted_before_boundary_filter & boundary_hit_mask);
localized.method_info.num_input_events = num_events;
localized.method_info.num_fitted_events = numel(planeset_subset);
end
