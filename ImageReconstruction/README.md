# ImageReconstruction

MATLAB UI and reconstruction pipeline for gamma-detector image reconstruction.

## Layout

- `ImageReconApp.m` - launch the UI.
- `IMAGE_default_config.m` - default project-local configuration.
- `IMAGE_run_recon.m` - reconstruction backend called by the UI.
- `method/` - localization methods used by the backend.
- `method/fabbri_analytical_lsf/` - Fabbri analytical fitting helpers.
- `library-construction/` - LSE reference libraries.
- `CALIB/` - optional CDF/UCM correction files.
- `input_data/` - default input folder for `.mat` files containing `planeset`.
- `IMAGE_Recon_Results/` - generated at runtime for reconstruction outputs.

## Start

```matlab
addpath('D:\110work\LineData-2022\gamma-detector\ImageReconstruction')
ImageReconApp
```

The UI can also process manually selected `.mat` files from any folder. Input
files must contain a `planeset` variable.
