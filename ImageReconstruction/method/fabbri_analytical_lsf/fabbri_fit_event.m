function fit = fabbri_fit_event(q, X, Y, cfg)
%FABBRI_FIT_EVENT Fit one 8 x 8 event to the analytical light spread model.

sum_q = sum(q(:));
if sum_q <= 0 || ~isfinite(sum_q)
    fit = empty_fit();
    return;
end

x_init = sum(q(:) .* X(:)) / sum_q;
y_init = sum(q(:) .* Y(:)) / sum_q;
amp_init = max(q(:));
if amp_init <= 0 || ~isfinite(amp_init)
    amp_init = sum_q / numel(q);
end

theta = [x_init; y_init; log(max(amp_init, 1e-12)); log(cfg.initial_alpha)];
theta = clamp_theta(theta, cfg, q);
lambda = cfg.initial_lambda;
converged = false;
accepted_any = false;

[residual, jacobian, obj] = model_residual(theta, q, X, Y, cfg);
last_delta = zeros(4, 1);

for iter = 1:cfg.max_iterations
    gradient = jacobian' * residual;
    hessian = jacobian' * jacobian;
    diag_h = abs(diag(hessian));
    diag_h(diag_h == 0) = 1;

    accepted = false;
    for trial = 1:cfg.max_inner_trials
        step_matrix = hessian + lambda * diag(diag_h);
        delta = -step_matrix \ gradient;

        if any(~isfinite(delta))
            lambda = lambda * 10;
            continue;
        end

        delta = limit_step(delta, cfg);
        theta_trial = clamp_theta(theta + delta, cfg, q);
        [res_trial, jac_trial, obj_trial] = model_residual(theta_trial, q, X, Y, cfg);

        if obj_trial <= obj
            accepted = true;
            accepted_any = true;
            break;
        end
        lambda = lambda * 4;
    end

    if ~accepted
        break;
    end

    last_delta = theta_trial - theta;
    theta = theta_trial;
    residual = res_trial;
    jacobian = jac_trial;

    objective_drop = obj - obj_trial;
    obj = obj_trial;
    lambda = max(lambda * 0.5, 1e-12);

    small_position_step = norm(last_delta(1:2)) < cfg.position_tolerance_mm;
    small_objective_drop = abs(objective_drop) < cfg.objective_tolerance * max(1, obj);
    if small_position_step && small_objective_drop
        converged = true;
        break;
    end
end

if ~converged && accepted_any && norm(last_delta(1:2)) < cfg.position_tolerance_mm
    converged = true;
end

[~, ~, final_obj, model] = model_residual(theta, q, X, Y, cfg);
fit.x_mm = theta(1);
fit.y_mm = theta(2);
fit.amplitude = exp(theta(3));
fit.alpha = exp(theta(4));
fit.residual_sse = final_obj;
fit.residual_mse = final_obj / numel(q);
fit.iterations = iter;
fit.converged = converged;
fit.initial_xy = [x_init, y_init];
fit.model = model;
end

function fit = empty_fit()
fit.x_mm = NaN;
fit.y_mm = NaN;
fit.amplitude = NaN;
fit.alpha = NaN;
fit.residual_sse = NaN;
fit.residual_mse = NaN;
fit.iterations = 0;
fit.converged = false;
fit.initial_xy = [NaN, NaN];
fit.model = [];
end

function [residual, jacobian, obj, model] = model_residual(theta, q, X, Y, cfg)
x0 = theta(1);
y0 = theta(2);
amplitude = exp(theta(3));
alpha = exp(theta(4));
b = cfg.shape_b;

dx = X - x0;
dy = Y - y0;
d2 = dx.^2 + dy.^2;
u = 1 + alpha * d2;
shape = u .^ (-b);
model = amplitude * shape;
residual_img = model - q;
residual = residual_img(:);
obj = sum(residual .^ 2);

u_pow = u .^ (-b - 1);
j_x = 2 * amplitude * b * alpha .* dx .* u_pow;
j_y = 2 * amplitude * b * alpha .* dy .* u_pow;
j_log_amp = model;
j_log_alpha = -amplitude * b * alpha .* d2 .* u_pow;
jacobian = [j_x(:), j_y(:), j_log_amp(:), j_log_alpha(:)];
end

function theta = clamp_theta(theta, cfg, q)
pos_min = cfg.position_bounds_mm(1);
pos_max = cfg.position_bounds_mm(2);
theta(1) = min(max(theta(1), pos_min), pos_max);
theta(2) = min(max(theta(2), pos_min), pos_max);

amp_max = max(10 * max(q(:)), 10);
theta(3) = min(max(theta(3), log(1e-12)), log(amp_max));
theta(4) = min(max(theta(4), log(cfg.alpha_bounds(1))), log(cfg.alpha_bounds(2)));
end

function delta = limit_step(delta, cfg)
pos_step = norm(delta(1:2));
if pos_step > cfg.max_position_step_mm
    delta(1:2) = delta(1:2) / pos_step * cfg.max_position_step_mm;
end
delta(3) = min(max(delta(3), -cfg.max_log_step), cfg.max_log_step);
delta(4) = min(max(delta(4), -cfg.max_log_step), cfg.max_log_step);
end

