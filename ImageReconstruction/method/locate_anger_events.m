function localized = locate_anger_events(planeset, locator, ~)
%LOCATE_ANGER_EVENTS Position events with standard Anger or RTP weighting.

num_events = numel(planeset);
raw_x = NaN(num_events, 1);
raw_y = NaN(num_events, 1);
total_signal = NaN(num_events, 1);
valid_event = false(num_events, 1);

centers = locator.sensor_centers_mm;

for ii = 1:num_events
    img = double(planeset{ii});
    if ~isequal(size(img), [8, 8])
        if numel(img) ~= 64
            continue;
        end
        img = reshape(img, 8, 8);
    end

    img(~isfinite(img)) = 0;
    if locator.clip_negative
        img(img < 0) = 0;
    end

    total = sum(img(:));
    if total <= 0
        continue;
    end

    x_proj = sum(img, 1);
    y_proj = sum(img, 2)';

    switch locator.position_weight_mode
        case {'standard', 'anger'}
            denom_x = total;
            denom_y = total;
            x_weight = x_proj;
            y_weight = y_proj;
        case {'rtp', 'projection_squared', 'quadratic_projection'}
            x_weight = x_proj .^ locator.projection_power;
            y_weight = y_proj .^ locator.projection_power;
            denom_x = sum(x_weight);
            denom_y = sum(y_weight);
        otherwise
            error('Unknown Anger position_weight_mode: %s', locator.position_weight_mode);
    end

    if denom_x <= 0 || denom_y <= 0
        continue;
    end

    raw_x(ii) = sum(x_weight .* centers) / denom_x;
    raw_y(ii) = sum(y_weight .* centers) / denom_y;
    total_signal(ii) = total;
    valid_event(ii) = true;
end

localized.raw_x_mm = raw_x(valid_event);
localized.raw_y_mm = raw_y(valid_event);
localized.residual = NaN(nnz(valid_event), 1);
localized.valid_mask = true(nnz(valid_event), 1);
localized.method_info.total_signal = total_signal(valid_event);
localized.method_info.position_weight_mode = locator.position_weight_mode;
localized.method_info.projection_power = locator.projection_power;
end
