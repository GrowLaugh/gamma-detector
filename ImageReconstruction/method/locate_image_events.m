function localized = locate_image_events(planeset, locator, cfg)
%LOCATE_IMAGE_EVENTS Dispatch planeset positioning to the selected method.

method = lower(strtrim(char(cfg.localization.method)));

switch method
    case {'lse', 'lse_softmax'}
        localized = locate_lse_softmax_events(planeset, locator, cfg);
    case {'lse64', 'lse_64ch', 'lse_softmax_64ch'}
        localized = locate_lse_softmax_64ch_events(planeset, locator, cfg);
    case {'anger', 'anger_standard', 'anger_rtp'}
        localized = locate_anger_events(planeset, locator, cfg);
    case {'fabbri', 'analytical_lsf', 'fabbri_analytical_lsf'}
        localized = locate_fabbri_analytical_lsf_events(planeset, locator, cfg);
    otherwise
        error('Unsupported localization method: %s', cfg.localization.method);
end
end
