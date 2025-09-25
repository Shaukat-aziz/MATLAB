# PPLN_simulation — Starter project for SPDC in PPLN

This project contains a compact MATLAB starter implementation for simulating SPDC in periodically-poled lithium niobate (PPLN), computing an interferogram R(Δl_i) and producing an animation that visualizes the optical layout.

## Files
- `params.m` — central parameters (edit to change physical values).
- `refractive_index.m` — Sellmeier-based refractive index (placeholder).
- `phase_match.m` — Δk and phase-matching amplitude.
- `two_photon_amplitude.m` — approximate JSI / weights along energy-conserving slice.
- `integrate_R0.m` — performs numeric integration and returns R(Δl_i).
- `sample_model.m` — sample transmission / phase model for idler-sample interaction.
- `simulation_driver.m` — top-level script: run this to compute and animate.
- `animate_beams.m` — creates an MP4 animation and a 2D cartoon.
- `plot_results.m` — plot utilities for analyzing results.
- `helper_math.m` — small helper functions.
- `export_video.m` — wrapper for video export.
- `outputs/` — directory where results are written.
- `params/` — optional folder for custom parameter `.mat` files.

## How to run (quick)
1. Open MATLAB and navigate to the project folder.
2. Run:
   ```matlab
   simulation_driver
