function locator = prepare_anger_locator(cfg)
%PREPARE_ANGER_LOCATOR Configure standard Anger / RTP positioning.

if ~isfield(cfg, 'anger')
    cfg.anger = struct();
end
if ~isfield(cfg.anger, 'sensor_centers_mm') || isempty(cfg.anger.sensor_centers_mm)
    cfg.anger.sensor_centers_mm = [-21, -15, -9, -3, 3, 9, 15, 21];
end
if ~isfield(cfg.anger, 'position_weight_mode') || isempty(cfg.anger.position_weight_mode)
    cfg.anger.position_weight_mode = 'rtp';
end
if ~isfield(cfg.anger, 'projection_power') || isempty(cfg.anger.projection_power)
    cfg.anger.projection_power = 2;
end
if ~isfield(cfg.anger, 'clip_negative') || isempty(cfg.anger.clip_negative)
    cfg.anger.clip_negative = true;
end

locator.method = 'anger';
locator.sensor_centers_mm = double(cfg.anger.sensor_centers_mm(:))';
method_name = lower(strtrim(char(cfg.localization.method)));
switch method_name
    case 'anger_standard'
        locator.position_weight_mode = 'standard';
    case 'anger_rtp'
        locator.position_weight_mode = 'rtp';
    otherwise
        locator.position_weight_mode = lower(strtrim(char(cfg.anger.position_weight_mode)));
end
locator.projection_power = cfg.anger.projection_power;
locator.clip_negative = cfg.anger.clip_negative;

switch locator.position_weight_mode
    case {'standard', 'anger'}
        locator.algorithm_name = 'Anger_Standard';
    case {'rtp', 'projection_squared', 'quadratic_projection'}
        locator.algorithm_name = sprintf('Anger_RTP_p%.3g', locator.projection_power);
    otherwise
        error('Unknown Anger position_weight_mode: %s', cfg.anger.position_weight_mode);
end
end
