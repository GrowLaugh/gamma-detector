function [coords, fitinfo] = fabbri_position_events(planeset, cfg)
%FABBRI_POSITION_EVENTS Fit all events in a planeset cell array.

[X, Y] = fabbri_anode_grid(cfg);
n_events = numel(planeset);

coords = NaN(n_events, 2, 'single');
residual_sse = NaN(n_events, 1, 'single');
residual_mse = NaN(n_events, 1, 'single');
iterations = zeros(n_events, 1, 'uint16');
converged = false(n_events, 1);
alpha = NaN(n_events, 1, 'single');
amplitude = NaN(n_events, 1, 'single');
initial_xy = NaN(n_events, 2, 'single');

for i = 1:n_events
    [q, ok] = fabbri_preprocess_event(planeset{i}, cfg);
    if ok
        fit = fabbri_fit_event(q, X, Y, cfg);
        coords(i, :) = single([fit.x_mm, fit.y_mm]);
        residual_sse(i) = single(fit.residual_sse);
        residual_mse(i) = single(fit.residual_mse);
        iterations(i) = uint16(fit.iterations);
        converged(i) = fit.converged;
        alpha(i) = single(fit.alpha);
        amplitude(i) = single(fit.amplitude);
        initial_xy(i, :) = single(fit.initial_xy);
    end

    if cfg.progress_interval > 0 && mod(i, cfg.progress_interval) == 0
        fprintf('    %d/%d events fitted\n', i, n_events);
    end
end

fitinfo.residual_sse = residual_sse;
fitinfo.residual_mse = residual_mse;
fitinfo.iterations = iterations;
fitinfo.converged = converged;
fitinfo.alpha = alpha;
fitinfo.amplitude = amplitude;
fitinfo.initial_xy = initial_xy;
end

