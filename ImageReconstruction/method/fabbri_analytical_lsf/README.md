# Analytical Fabbri least-square fitting

Minimal Fabbri analytical fitting functions used by the IMAGE production
pipeline.

Call chain:

1. `locate_fabbri_analytical_lsf_events.m`
2. `fabbri_position_events.m`
3. `fabbri_preprocess_event.m`
4. `fabbri_fit_event.m`

Supporting helpers:

- `fabbri_default_config.m`
- `fabbri_anode_grid.m`

The historical top-level `Fabbri/` folder contains older diagnostics and
standalone plotting scripts. Those files are not required for
`IMAGE_main_recon.m`.

