% Calculate daily ASI regimes from OAFlux turbulent heat flux and SST,
% 1985-2020 (main results, Figs. 2-8).
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  P.oaflux_NIO      (pipelines/asi_regimes_oaflux_processing.m)
% Outputs: P.corrs_OAFlux    lagged SST-Qturb and dSST/dt-Qturb correlations
%          P.asi_OAFlux      ASI regime parameter P_o (asi_regimes.param)

clear; clc

P = asi_regimes_paths;

%% Load OAFlux

load(P.oaflux_NIO, 'flux_turb_NIO')

dt = datetime(datevec(flux_turb_NIO.time));
[~, keep_noleap] = doy_noleap(dt);

% OAFlux latent and sensible heat fluxes are positive upward
Qturb = flux_turb_NIO.lhf(:,:,keep_noleap) + flux_turb_NIO.shf(:,:,keep_noleap);
SST = flux_turb_NIO.sst(:,:,keep_noleap);

SST(SST > 50) = nan;
Qturb(isnan(SST)) = nan;

%% Calculate ASI regimes

[sst_thf_anomaly_corrs, asi_regimes] = asi_regimes_calculation( ...
    SST, Qturb, dt(keep_noleap), flux_turb_NIO.lon, flux_turb_NIO.lat);

sst_thf_anomaly_corrs.source = 'OAFlux';
asi_regimes.source = 'OAFlux';

%% Save

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.corrs_OAFlux, 'sst_thf_anomaly_corrs', '-v7.3')
save(P.asi_OAFlux, 'asi_regimes', '-v7.3')
