function locator = prepare_lse_softmax_64ch_locator(cfg)
%PREPARE_LSE_SOFTMAX_64CH_LOCATOR Load the full 64-channel LSE library on GPU.

if parallel.gpu.GPUDevice.isAvailable
    g = gpuDevice;
    fprintf('Using GPU for 64-channel LSE localization: %s\n', g.Name);
    if isfield(cfg.localization, 'reset_gpu') && cfg.localization.reset_gpu
        reset(g);
    end
else
    error('No GPU detected. 64-channel LSE softmax localization uses gpuArray.');
end

if ~exist(cfg.library_file, 'file')
    error('LSE library file not found: %s', cfg.library_file);
end

fprintf('Loading 64-channel LSE reference library...\n');
lib_data = load(cfg.library_file, 'lightMapLibrary');
if ~isfield(lib_data, 'lightMapLibrary')
    error('64-channel LSE localization requires lightMapLibrary in: %s', cfg.library_file);
end

raw_lib_cell = lib_data.lightMapLibrary;
[Lx, Ly, ~] = size(raw_lib_cell);
dx_lib = (cfg.lib_phys_end - cfg.lib_phys_start) / (Lx - 1);
library_tensor = zeros(Lx, Ly, 64, 'single');
use_softmax_interpolation = true;
if isfield(cfg.localization, 'enable_softmax_interpolation')
    use_softmax_interpolation = logical(cfg.localization.enable_softmax_interpolation);
end

for ix = 1:Lx
    for iy = 1:Ly
        if ndims(raw_lib_cell) <= 2
            response = raw_lib_cell{ix, iy};
        else
            response = raw_lib_cell{ix, iy, 1};
        end

        if isempty(response) || numel(response) ~= 64
            vec64 = ones(1, 64, 'single') * 1e-6;
        else
            vec64 = reshape(single(response), 1, 64);
            vec64(~isfinite(vec64)) = 0;
            vec64(vec64 < 0) = 0;
        end

        library_tensor(ix, iy, :) = reshape(vec64 ./ (sum(vec64) + eps('single')), 1, 1, 64);
    end
end

locator.method = 'lse_softmax_64ch';
if use_softmax_interpolation
    locator.algorithm_name = 'LSE_Gaussian_64ch';
else
    locator.algorithm_name = 'LSE_Top1_64ch';
end
locator.use_softmax_interpolation = use_softmax_interpolation;
locator.Lx = Lx;
locator.Ly = Ly;
locator.dx_lib = dx_lib;
locator.G_LibTensor64 = gpuArray(library_tensor);
end
