# asi_regimes

MATLAB code accompanying

> Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air–Sea Interaction Reveals Ocean Control of Turbulent Heat Flux in the South Asian Summer Monsoon. *J. Climate*, **39**, 6319–6333, https://doi.org/10.1175/JCLI-D-26-0118.1.

The code runs the idealized stochastic model, processes the observational and reanalysis data, computes the air–sea interaction (ASI) regime parameter *P<sub>o</sub>*, and generates the figures, supplemental figure and supplemental video. Please cite the paper if you use this code.

## Repository layout

```
asi_regimes_paths.m   all data, intermediate-file and figure paths (edit this first)
pipelines/            data processing, the P_o calculation, the Monte Carlo test and the stochastic model
figures/              one script per figure (Figs. 1–8, Fig. S1, Video 1)
analysis/             numbers quoted in the text (Sections 3a and 4b)
util/                 shared functions and bundled third-party code
```

Each script starts with a header that lists its inputs and outputs. Every script calls `P = asi_regimes_paths;`, which returns the file locations and adds `util/` and `util/gsw/` to the MATLAB path, so scripts can be run from any folder.

## Requirements

- MATLAB R2023a or later (`xregion` in Fig. 8). The code was tested with R2025b and R2026a.
- Statistics and Machine Learning Toolbox (`corr`, `prctile`, `tcdf`).
- Mapping Toolbox (`shaperead`, for the coastlines in Figs. 3 and 5).
- Memory: the OAFlux pipelines run comfortably in 16 GB. The ERA5 parameter calculation (0.25° grid) needs tens of GB.

No other toolboxes or File Exchange downloads are needed; the third-party code that is used is bundled in `util/` (see below).

## Data

Download the following datasets, then set the four user paths at the top of `asi_regimes_paths.m`:

| Dataset | Source | Path setting | Expected files |
|---|---|---|---|
| OAFlux daily turbulent fluxes (with the OISST v1 SST distributed with OAFlux), 1985–2020 | Asia-Pacific Data Research Center, https://apdrc.soest.hawaii.edu/datadoc/whoi_oaflux.php | `P.oaflux_dir` | `turbulence/{lh,sh,ts,ta,qa,ws}_oaflux_<year>.nc`, `evaporation/evapr_oaflux_<year>.nc` |
| ERA5 reanalysis surface sensible and latent heat flux and SST, 1985–2020, Indian Ocean | Copernicus Climate Data Store, https://cds.climate.copernicus.eu/datasets/reanalysis-era5-complete | `P.era5_dir` | `ERA5_THF_SST_IO_<year>_rates.nc`, 6-hourly, with `msshf`, `mslhf` and `sst` |
| CERES SYN1deg-Day Ed4.1 surface fluxes, March 2000 – March 2024, 40–100°E, 10°S–30°N | NASA Langley, https://ceres.larc.nasa.gov/data/ | `P.data_dir` | `CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20000301-20240331.nc`, with the adjusted all-sky surface SW and LW up and down fluxes |
| Arabian Sea Experiment mooring (1994–95) | WHOI Upper Ocean Processes Group, https://uop.whoi.edu/projects/arabiansea/arabiansea.html | `P.data_dir` | `arabianseaflux_1hr.epic`, `arabianseatemp.nc.epic`, `arabianseasal.nc.epic`, `arabianseavel.epic` |
| ASIRI Bay of Bengal mooring | WHOI Upper Ocean Processes Group, https://uop.whoi.edu/projects/Bengal/QCData.html | `P.data_dir` | `beng1_flux_1hr.mat`, `bob1_zsal_hr_vfilled.mat`, `bob1_zcur_hr.mat` |

`P.output_dir` is where intermediate `.mat` files are written and read. `P.fig_dir` (by default `output_figures/` in the repository) is where figures and the video are written.

## Run order

Run the pipelines first; each writes a `.mat` file to `P.output_dir` that later scripts read. Run times are approximate, measured on a 16-core, 48 GB workstation.

| Step | Script | Produces | Approx. time |
|---|---|---|---|
| 1 | `pipelines/asi_regimes_oaflux_processing.m` | `oaflux_NIO.mat` | minutes |
| 2 | `pipelines/asi_regimes_era5_processing.m` | `ERA5_THF_SST_IO_1985_2020.mat` | not timed |
| 3 | `pipelines/asi_regimes_ceres_processing.m` | `CERES_rad.mat` | minutes |
| 4 | `pipelines/asi_regimes_parameter_calculation.m` | `sst_thf_anomaly_corrs_OAFlux.mat`, `asi_regimes_OAFlux.mat` | ~1.5–2 h |
| 5 | `pipelines/asi_regimes_parameter_calculation_ERA5.m` | `sst_thf_anomaly_corrs_ERA5.mat`, `asi_regimes_ERA5.mat` | several hours (not timed) |
| 6 | `pipelines/monte_carlo_resampling.m` | `sst_thf_anomaly_corrs_monte_carlo.mat` | ~1 h |
| 7 | `pipelines/FCL98_model.m` | `FCL_runs_tropical_params.mat` | ~1.5 h |

Steps 2, 3 and 5 are needed only for Fig. 6, Video 1 and the radiative-flux analysis. Step 7 is needed only for Figs. 1 and 8.

Then run the figure and analysis scripts in any order:

| Paper item | Script | Needs steps |
|---|---|---|
| Fig. 1 | `figures/asi_regimes_fig1_curves.m` | 7 |
| Fig. 2 | `figures/asi_regimes_fig2_decomposition.m` | 1, 4 |
| Fig. 3 | `figures/asi_regimes_fig3_seasonal.m` | 1, 4 |
| Fig. 4 and the percentages in Section 3a | `figures/asi_regimes_fig4_histograms.m` | 1, 4 |
| Fig. 5 | `figures/asi_regimes_fig5_wind.m` | 1, 4 |
| Fig. 6, the OAFlux–ERA5 correlation (Section 3c) and Video 1 | `figures/asi_regimes_fig6_oaflux_vs_era5.m` (set `make_video = true` for the video) | 1, 2, 4, 5 |
| Fig. 7 and Fig. S1 | `figures/asi_regimes_fig7_mooring.m` | 1, 4 |
| Fig. 8 | `figures/asi_regimes_fig8_simple_model.m` | 1, 4, 7 |
| Nearshore vs. offshore *P<sub>o</sub>* (Section 3a) | `analysis/asi_regimes_nearshore_vs_offshore.m` | 1, 4 |
| SST tendency–radiative flux correlations (Section 4b) | `analysis/asi_regimes_radiative_correlation.m` | 1, 3 |
| Monte Carlo noise threshold (appendix) | printed by `pipelines/monte_carlo_resampling.m` | 1, 4 |

## Method summary

`util/asi_regimes_calculation.m` implements Section 2b and is shared by the OAFlux and ERA5 pipelines:

1. Anomalies (`util/asi_anomalies_noleap.m`): data are put on a 365-day no-leap calendar. The daily climatology uses an 11-day circular day-of-year window and is smoothed with a 31-day circular running mean. SST tendency is computed from the SST anomalies.
2. Lagged correlations (`util/lagged_corrs_doy.m`): SST–Q<sub>turb</sub> and dSST/dt–Q<sub>turb</sub> correlations at lags of −30 to 30 days, for each day of the year, from all years within a 31-day circular window.
3. Decomposition: each pair of curves is split into symmetric and antisymmetric parts and fitted with a weighted nonnegative least-squares fit, which gives *P<sub>o</sub>*. The weights use effective-sample-size p values, Benjamini–Hochberg q values (`util/bh_qvalues.m`) and a sigmoid.

## Differences from the published paper

The code was cleaned up for release, and a few results differ slightly from the paper. None changes the paper's conclusions.

- **Monte Carlo test (appendix).** The released script uses the same no-leap, circular-window method as the main calculation, with a fixed random seed. It gives a 95th-percentile threshold of 0.186, which flags 0.15% of the observed curves as noise. The paper reports 0.194 and 0.23%. Either way, very few curves are consistent with noise.
- **Stochastic model (Figs. 1 and 8).** The model is now seeded (`rng(0)`), so reruns are reproducible. The published figures came from an earlier, unseeded run whose forcing set was ω<sub>a</sub> = [0, 0.15, linspace(0.3, 10, 60)]. The released script uses ω<sub>a</sub> = `linspace(0, 10, 60)`, as in Table 1. Fig. 1 is visually unchanged. The *P<sub>o</sub>* curve in Fig. 8 differs slightly because the noise realization is different.
- **Fig. 8 shaded band.** The band is now computed in the script as the 5th–95th percentile range of the monthly climatological *P<sub>o</sub>* over the whole domain. This gives *P<sub>o</sub>* = 0.123–0.976. The published figure used 0.138–0.949, values carried over from an earlier draft that cannot be reproduced. The right edge of the band therefore moves from a forcing ratio of about 2.6 to about 3.
- **Fig. 4.** The histograms and percentages use all grid points in the domain, as in the published figure and text. The published caption says that only the Arabian Sea and Bay of Bengal are included.
- **Fig. 5.** The script reproduces the figure. The highest correlation it prints is R = 0.54; the text reports 0.55.
- **Fig. 1.** The script produces the panel plots. The annotations in the published figure were added in a graphics editor.
- **Figure layout.** Figures request the published window sizes in pixels, but font scaling and window limits depend on the screen and MATLAB version, so sizes and spacing can differ from the published figures. The data shown are the same.
- **Radiative-flux correlations (Section 4b).** The analysis now uses the same anomaly method as the main calculation, and the OAFlux SST over the full CERES overlap (March 2000 – December 2020). Of the significant simultaneous correlations, the median is 0.18 and the distribution peaks near 0.15, below 0.2, and drops off steeply above it, as stated in the paper.

## Third-party code in `util/`

- `cmocean.m`: cmocean colormaps (Thyng et al. 2016, *Oceanography*, https://doi.org/10.5670/oceanog.2016.66), MATLAB version by Chad Greene. MIT license, in `util/cmocean_LICENSE`.
- `gsw/`: the subset of the Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10, v3.06.13, needed to compute seawater density in Fig. 7 (McDougall and Barker 2011, https://www.teos-10.org). License in `util/gsw/gsw_licence.html`.
- `ne_10m_coastline/`, `ne_50m_coastline/`: Natural Earth coastlines (https://www.naturalearthdata.com), public domain.

## License

MIT; see `LICENSE`. The bundled third-party code keeps its own license, as listed above.
