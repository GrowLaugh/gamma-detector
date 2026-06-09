function [X, Y] = fabbri_anode_grid(cfg)
%FABBRI_ANODE_GRID Return physical coordinates for the 8 x 8 anode centers.

centers = ((1:8) - 4.5) * cfg.anode_pitch_mm;
[X, Y] = meshgrid(centers, centers);
end

