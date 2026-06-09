# gamma-detector

Gamma detector data processing and image reconstruction toolkit.

This repository contains three main modules:

1. Monte Carlo simulation library construction: `libraryGenerate/`
2. Raw detector data preprocessing: `rawDataProcess/`
3. Image reconstruction: `ImageReconstruction/` **(main module)**

The typical workflow is:

```text
Monte Carlo hits / reference data
        -> libraryGenerate
        -> LSE reference library

Raw detector .bin files
        -> rawDataProcess
        -> calibrated event .mat files with planeset
        -> ImageReconstruction
        -> reconstructed images, tables, and analysis outputs
```

## Repository Layout

```text
gamma-detector/
+-- libraryGenerate/        # Monte Carlo reference-library construction
+-- rawDataProcess/         # Raw .bin preprocessing and ADC calibration
+-- ImageReconstruction/    # MATLAB UI and reconstruction backend
+-- MLE/                    # MLE-related legacy/experimental code
+-- README.md
```

## 1. Monte Carlo Library Construction

Directory:

```text
libraryGenerate/
```

This module builds detector response libraries from Monte Carlo simulation output. The generated library is used by LSE-based image reconstruction methods.

Main files:

```text
parameters.json
newbuildreferencelibrary.m
build_2d_library_from_3d.m
build_symmetric_library.m
```

`parameters.json` defines the simulation grid, detector geometry, input hit-data folder, and output library name. Example fields include:

```json
{
  "grid_parameters": {
    "x": {"start": 2.0, "end": 23.0, "step": 3},
    "y": {"start": 2.0, "end": 23.0, "step": 3},
    "z": {"start": -6.5, "end": -0.5, "step": 1.0}
  },
  "detector_parameters": {
    "pixels_x": 8,
    "pixels_y": 8,
    "pixel_size_mm": 6.0625
  }
}
```

Typical use:

```matlab
cd libraryGenerate
newbuildreferencelibrary
```

Expected output is a `.mat` reference library containing variables such as:

```text
lightMapLibrary
referprojection
```

For reconstruction, place the generated `.mat` library under:

```text
ImageReconstruction/library-construction/
```

## 2. Raw Data Preprocessing

Directory:

```text
rawDataProcess/
```

This module converts raw detector binary data into calibrated event `.mat` files. These `.mat` files are the input of the image reconstruction module.

Main files:

```text
correction.m              # Generate Offsetdata.mat and Gaindata.mat
step1_mouse_all.m         # Batch process .bin files into calibrated events
step2_anger_multiple.m    # Anger-positioning legacy batch workflow
Disp.m                    # Visualization helper
Offsetdata.mat            # Example/default offset calibration
Gaindata.mat              # Example/default gain calibration
```

### ADC Calibration

Run `correction.m` after editing the raw data folder and file names:

```matlab
cd rawDataProcess
correction

```

It reads offset/gain `.bin` files and produces:

```text
Offsetdata.mat
Gaindata.mat

```

### Event Generation

Run `step1_mouse_all.m` after editing:

```matlab
input_data_folder
file_basename
file_indices_to_process
output_folder
energy_spec_folder

```

The output event files should contain a variable named:

```matlab
planeset

```

Those files can then be copied or selected directly in the reconstruction UI.

## 3. Image Reconstruction

Directory:

```text
ImageReconstruction/

```

This is the main module of the repository. It provides a MATLAB UI for selecting input files, localization methods, post-processing mode, reconstruction parameters, and output location.

### Key Files

```text
ImageReconApp.m           # MATLAB UI entry point
IMAGE_default_config.m    # Default project-local configuration
IMAGE_run_recon.m         # Reconstruction backend called by the UI
method/                   # Localization algorithms
library-construction/     # LSE reference libraries
CALIB/                    # Optional CDF/UCM correction files
input_data/               # Default input folder for event .mat files

```

### Start the UI

From MATLAB:

```matlab
addpath('D:\110work\LineData-2022\gamma-detector\ImageReconstruction')
ImageReconApp

```

Or, after changing to the project folder:

```matlab
cd('D:\110work\LineData-2022\gamma-detector\ImageReconstruction')
ImageReconApp

```

### Input Data

The reconstruction backend expects `.mat` files containing:

```matlab
planeset

```

You can either:

- place files under `ImageReconstruction/input_data/`, then use pattern + indices mode;
- or use manual file selection in the UI, which is recommended when file names are irregular.

When using manual file selection, output names preserve the input file name. For example:

```text
array10by10lead060503.mat

```

may produce:

```text
Recon_LSE_Gaussian_64ch_gridmask_array10by10lead060503.png
ReconData_LSE_Gaussian_64ch_gridmask_array10by10lead060503.mat

```

### Localization Methods

The UI currently supports:

```text
lse_softmax_64ch
lse_softmax
fabbri_analytical_lsf
anger_standard
anger_RTP

```

LSE methods require a reference library in:

```text
ImageReconstruction/library-construction/

```

Fabbri analytical fitting uses helper functions in:

```text
ImageReconstruction/method/fabbri_analytical_lsf/

```

### Post-processing Modes

The UI supports:

```text
flood
gridmask
slit
single_hole
single_hole_scan

```

Mode summary:

- `flood`: reconstructs flood image and saves X/Y projections.
- `gridmask`: reconstructs grid-mask data and extracts peak/valley profile information.
- `slit`: builds aligned slit profiles and estimates FWHM.
- `single_hole`: normalizes and accumulates single-hole images.
- `single_hole_scan`: analyzes position linearity across a scan sequence.

### Optional Corrections

CDF and UCM correction files are included under:

```text
ImageReconstruction/CALIB/Calibration_Results_FloodCDF_LSE/
ImageReconstruction/CALIB/Calibration_Results_AngerCDF/

```

They are only used when the corresponding UI checkboxes are enabled.

### Output

By default, reconstruction outputs are generated under:

```text
ImageReconstruction/IMAGE_Recon_Results/

```

Typical subfolders:

```text
reconstruction/       # main reconstruction PNG/MAT files
flood/                # flood projections
gridmask/             # gridmask profiles and tables
slit/                 # aligned slit profiles and summaries
single_hole/          # per-file and composite images
single_hole_scan/     # linearity figures and tables

```

## MATLAB Requirements

Recommended environment:

- MATLAB with App Designer/uifigure support
- Image Processing / plotting support used by standard MATLAB graphics
- Parallel Computing Toolbox and CUDA-capable GPU for LSE methods using `gpuArray`

If no GPU is available, use `anger_standard`, `anger_RTP`, or `fabbri_analytical_lsf` instead of LSE methods.

## Notes

- `rawDataProcess` scripts currently contain editable path variables near the top of each script. Update those paths before running preprocessing.
- `ImageReconstruction/input_data/` is intentionally kept mostly empty. Put generated event `.mat` files there or select them manually in the UI.
- `ImageReconstruction/IMAGE_Recon_Results/` is runtime output and should not be committed unless specific result examples are needed.
