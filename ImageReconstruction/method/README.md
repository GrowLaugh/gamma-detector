# IMAGE method folder

This folder contains event-localization methods used by `IMAGE_main_recon.m`.
The main script handles reconstruction, rendering, and post-processing; this
folder should only contain method setup and event-positioning code.

## Dispatchers

- `prepare_image_locator.m`: selects and initializes the requested method.
- `locate_image_events.m`: dispatches `planeset` events to the selected method.

## Anger centroid / RTP

- `prepare_anger_locator.m`: configures sensor centers and Anger/RTP mode.
- `locate_anger_events.m`: converts each event to an Anger-position estimate.

## LSE softmax

- `prepare_lse_softmax_locator.m`: loads the MC light-map library and prepares
  the GPU reference tensor.
- `locate_lse_softmax_events.m`: runs local softmax-weighted LSE localization.
  Set `cfg.localization.enable_softmax_interpolation = false` to disable
  top-k interpolation and use the top-1 minimum-LSE library point directly.

## LSE softmax 64-channel

- `prepare_lse_softmax_64ch_locator.m`: loads the full 8 x 8 MC light-map
  library and prepares the GPU reference tensor with 64 anode channels.
- `locate_lse_softmax_64ch_events.m`: runs the same coarse search and local
  softmax-weighted LSE localization using all 64 normalized anode values.
  Select it with `cfg.localization.method = 'lse_softmax_64ch'`.

## Analytical Fabbri least-square fitting

- `prepare_fabbri_analytical_lsf_locator.m`: loads the local Fabbri dependency
  folder and prepares fitting parameters.
- `locate_fabbri_analytical_lsf_events.m`: runs analytical fitting for a
  `planeset`.
- `fabbri_analytical_lsf/`: minimal Fabbri fitting implementation required by
  the production IMAGE workflow.
