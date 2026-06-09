function locator = prepare_lse_softmax_locator(cfg)
%PREPARE_LSE_SOFTMAX_LOCATOR Load the LSE reference library on GPU.

if parallel.gpu.GPUDevice.isAvailable
    g = gpuDevice;
    fprintf('Using GPU for LSE localization: %s\n', g.Name);
    if isfield(cfg.localization, 'reset_gpu') && cfg.localization.reset_gpu
        reset(g);
    end
else
    error('No GPU detected. LSE softmax localization uses gpuArray.');
end

if ~exist(cfg.library_file, 'file')
    error('LSE library file not found: %s', cfg.library_file);
end

fprintf('Loading LSE reference library...\n');
lib_data = load(cfg.library_file, 'lightMapLibrary');
if isfield(lib_data, 'lightMapLibrary')
    raw_lib_cell = lib_data.lightMapLibrary;
else
    lib_data = load(cfg.library_file, 'referprojection');
    raw_lib_cell = lib_data.referprojection;
end

[Lx, Ly, ~] = size(raw_lib_cell);
dx_lib = (cfg.lib_phys_end - cfg.lib_phys_start) / (Lx - 1);
library_tensor = zeros(Lx, Ly, 16, 'single');
use_softmax_interpolation = true;
if isfield(cfg.localization, 'enable_softmax_interpolation')
    use_softmax_interpolation = logical(cfg.localization.enable_softmax_interpolation);
end

for ix = 1:Lx
    for iy = 1:Ly
        if ndims(raw_lib_cell) <= 2
            tmp = raw_lib_cell{ix, iy};
        else
            tmp = raw_lib_cell{ix, iy, 1};
        end

        if isempty(tmp)
            vec16 = ones(1, 1, 16, 'single') * 1e-6;
        else
            if numel(tmp) == 64
                img = single(tmp);
                s = sum(img(:));
                if s > 0
                    img = img / s;
                end
                vec16 = reshape([sum(img, 1), sum(img, 2)'], 1, 1, 16);
            elseif numel(tmp) == 16
                vec = single(tmp(:))';
                s = sum(vec);
                if s > 0
                    vec = vec / s;
                end
                vec16 = reshape(vec, 1, 1, 16);
            else
                vec16 = ones(1, 1, 16, 'single') * 1e-6;
            end
        end

        library_tensor(ix, iy, :) = vec16 / (sum(vec16(:)) + eps);
    end
end

locator.method = 'lse_softmax';
if use_softmax_interpolation
    locator.algorithm_name = 'LSE_Gaussian';
else
    locator.algorithm_name = 'LSE_Top1';
end
locator.use_softmax_interpolation = use_softmax_interpolation;
locator.Lx = Lx;
locator.Ly = Ly;
locator.dx_lib = dx_lib;
locator.G_LibTensor = gpuArray(library_tensor);
end
