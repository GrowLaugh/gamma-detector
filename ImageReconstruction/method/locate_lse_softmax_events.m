function localized = locate_lse_softmax_events(planeset, locator, cfg)
%LOCATE_LSE_SOFTMAX_EVENTS Position events with the current LSE softmax method.

num_events = numel(planeset);
all_evt_16 = zeros(num_events, 16, 'single');

for ii = 1:num_events
    all_evt_16(ii, :) = extract_feature_16ch(planeset{ii});
end

final_x_mm = zeros(num_events, 1, 'single');
final_y_mm = zeros(num_events, 1, 'single');
residual = zeros(num_events, 1, 'single');
if isfield(locator, 'use_softmax_interpolation')
    use_softmax_interpolation = locator.use_softmax_interpolation;
else
    use_softmax_interpolation = true;
end

num_batches = ceil(num_events / cfg.localization.batch_size);
for bb = 1:num_batches
    curr_idx = (bb-1) * cfg.localization.batch_size + 1 : ...
        min(bb * cfg.localization.batch_size, num_events);
    batch_evt = all_evt_16(curr_idx, :);

    [rx, ry] = batch_rough_scheme2(batch_evt, locator.Lx, locator.Ly, ...
        cfg.detector_pitch_mm, cfg.lib_phys_start, locator.dx_lib);

    G_evt = gpuArray(batch_evt);
    G_rx = gpuArray(single(rx));
    G_ry = gpuArray(single(ry));

    [G_fx, G_fy, G_min_lse] = gpu_vectorized_lse_baseline_opt(G_evt, ...
        locator.G_LibTensor, G_rx, G_ry, cfg.localization.coarse_radius, ...
        cfg.localization.lse_temperature, cfg.localization.baseline_ratio, ...
        cfg.localization.top_k_ratio, use_softmax_interpolation);

    final_x_mm(curr_idx) = (gather(G_fx) - 1) * locator.dx_lib + cfg.lib_phys_start;
    final_y_mm(curr_idx) = (gather(G_fy) - 1) * locator.dx_lib + cfg.lib_phys_start;
    residual(curr_idx) = gather(G_min_lse);
end

localized.raw_x_mm = double(final_x_mm);
localized.raw_y_mm = double(final_y_mm);
localized.residual = double(residual);
localized.valid_mask = true(num_events, 1);
end

function v16 = extract_feature_16ch(event_img)
img = single(event_img);
s_evt = sum(img(:));
if s_evt > 0
    img = img / s_evt;
end
xp = sum(img, 1);
yp = sum(img, 2)';
v16 = [xp, yp];
v16 = v16 / (sum(v16) + eps);
end

function [rx, ry] = batch_rough_scheme2(batch_evt, Lx, Ly, pixel_size_mm, lib_start, dx_lib)
xp = batch_evt(:, 1:8);
yp = batch_evt(:, 9:16);
cx = sum(xp.^2 .* (1:8), 2) ./ (sum(xp.^2, 2) + eps);
cy = sum(yp.^2 .* (1:8), 2) ./ (sum(yp.^2, 2) + eps);
x_phys = (cx - 4.5) * pixel_size_mm;
y_phys = (cy - 4.5) * pixel_size_mm;
rx = round((x_phys - lib_start) / dx_lib + 1);
ry = round((y_phys - lib_start) / dx_lib + 1);
rx = max(1, min(Lx, rx));
ry = max(1, min(Ly, ry));
end

function [fx, fy, min_lse] = gpu_vectorized_lse_baseline_opt(G_evt, G_lib, G_rx, G_ry, radius, temp, base_ratio, top_k_ratio, use_softmax_interpolation)
[batch_size, ~] = size(G_evt);
[Lx, Ly, ~] = size(G_lib);
[GX, GY] = meshgrid(-radius:radius, -radius:radius);
gx = gpuArray(single(GX(:)'));
gy = gpuArray(single(GY(:)'));

TX = bsxfun(@plus, G_rx, gx);
TY = bsxfun(@plus, G_ry, gy);
out_of_bounds_mask = (TX < 1) | (TX > Lx) | (TY < 1) | (TY > Ly);
TX_safe = max(1, min(Lx, TX));
TY_safe = max(1, min(Ly, TY));
G_idx = (TY_safe - 1) * Lx + TX_safe;

lib_flat = reshape(G_lib, [], 16);
ref_block = reshape(lib_flat(G_idx(:), :), batch_size, [], 16);
diff_block = bsxfun(@minus, ref_block, reshape(G_evt, batch_size, 1, 16));
lse = sum(diff_block.^2, 3);
lse_safe = lse + single(out_of_bounds_mask) * 1e6;

[min_lse, best_idx] = min(lse_safe, [], 2);
if ~use_softmax_interpolation
    dx = reshape(gx(best_idx), [], 1);
    dy = reshape(gy(best_idx), [], 1);
    fx = G_rx + dx;
    fy = G_ry + dy;
    return;
end

num_candidates = size(lse_safe, 2);
k_target = max(1, round(num_candidates * top_k_ratio));
k_next = min(k_target + 1, num_candidates);

sorted_lse = sort(lse_safe, 2, 'ascend');
threshold_lse = sorted_lse(:, k_target);
next_lse = sorted_lse(:, k_next);

w_k = exp(-(threshold_lse - min_lse) / temp);
w_k_next = exp(-(next_lse - min_lse) / temp);
condition_mask = single(w_k > 0.20);
adaptive_baseline = condition_mask .* w_k_next + (~condition_mask) .* base_ratio;

top_k_mask = lse_safe <= threshold_lse;
weights = exp(-(lse_safe - min_lse) / temp);
weights = weights .* single(top_k_mask);
weights = max(0, bsxfun(@minus, weights, adaptive_baseline));

sum_w = sum(weights, 2) + eps;
dx = sum(bsxfun(@times, weights, gx), 2) ./ sum_w;
dy = sum(bsxfun(@times, weights, gy), 2) ./ sum_w;

fx = G_rx + dx;
fy = G_ry + dy;
end
