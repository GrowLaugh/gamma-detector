function locator = prepare_image_locator(cfg)
%PREPARE_IMAGE_LOCATOR Build method-specific state for event positioning.

method = lower(strtrim(char(cfg.localization.method)));

switch method
    case {'lse', 'lse_softmax'}
        locator = prepare_lse_softmax_locator(cfg);
    case {'lse64', 'lse_64ch', 'lse_softmax_64ch'}
        locator = prepare_lse_softmax_64ch_locator(cfg);
    case {'anger', 'anger_standard', 'anger_rtp'}
        locator = prepare_anger_locator(cfg);
    case {'fabbri', 'analytical_lsf', 'fabbri_analytical_lsf'}
        locator = prepare_fabbri_analytical_lsf_locator(cfg);
    otherwise
        error('Unsupported localization method: %s', cfg.localization.method);
end
end
