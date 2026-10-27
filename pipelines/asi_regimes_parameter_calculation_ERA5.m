% Calculate daily ASI regimes from ERA5 turbulent heat flux and SST,
% 1985-2020 (dataset comparison, Fig. 6 and Video 1).
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  P.era5_1985_2020  (pipelines/asi_regimes_era5_processing.m)
% Outputs: P.corrs_ERA5      lagged SST-Qturb and dSST/dt-Qturb correlations
%          P.asi_ERA5        ASI regime parameter P_o (asi_regimes.param)
%
% The 0.25 degree ERA5 grid needs substantially more memory than OAFlux
% (tens of GB).

clear; clc

P = asi_regimes_paths;

%% Load ERA5

load(P.era5_1985_2020, 'ERA5_THF_SST_IO')

dt = datetime(ERA5_THF_SST_IO.dvec);
[~, keep_noleap] = doy_noleap(dt);

% ERA5 fluxes are positive downward; flip to positive upward
Qturb = -ERA5_THF_SST_IO.lhf(:,:,keep_noleap) + ...
        -ERA5_THF_SST_IO.shf(:,:,keep_noleap);
SST = ERA5_THF_SST_IO.sst(:,:,keep_noleap);

lon = ERA5_THF_SST_IO.lon;
lat = ERA5_THF_SST_IO.lat;
clear ERA5_THF_SST_IO

%% Calculate ASI regimes

[sst_thf_anomaly_corrs, asi_regimes] = asi_regimes_calculation( ...
    SST, Qturb, dt(keep_noleap), lon, lat);

sst_thf_anomaly_corrs.source = 'ERA5';
asi_regimes.source = 'ERA5';

%% Save

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.corrs_ERA5, 'sst_thf_anomaly_corrs', '-v7.3')
save(P.asi_ERA5, 'asi_regimes', '-v7.3')
