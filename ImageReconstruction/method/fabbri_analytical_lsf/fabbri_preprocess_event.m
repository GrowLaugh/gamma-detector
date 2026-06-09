function [q, ok] = fabbri_preprocess_event(raw_event, cfg)
%FABBRI_PREPROCESS_EVENT Convert an event to a finite 8 x 8 charge image.

q = double(raw_event);
if numel(q) ~= 64
    ok = false;
    q = zeros(8, 8);
    return;
end

if ~isequal(size(q), [8, 8])
    q = reshape(q, 8, 8);
end

q(~isfinite(q)) = 0;
if cfg.clip_negative
    q(q < 0) = 0;
end

total = sum(q(:));
ok = total > 0;
if ok && cfg.normalize_each_event
    q = q / total;
end
end

