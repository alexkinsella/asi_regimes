% Test whether the observed lagged correlation curves could arise from
% noise, using Monte Carlo resampling (appendix, Eq. A1).
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Synthetic SST, SST tendency, and Qturb anomaly series are drawn at
% random from the pooled OAFlux anomalies (all grid points and days),
% which destroys any real relationship between them. Lagged correlations
% are computed exactly as for the observations, and the 95th percentile
% of max|C_Ts,Q| + max|C_dTs/dt,Q| defines the noise threshold.
%
% Inputs:  P.oaflux_NIO    (pipelines/asi_regimes_oaflux_processing.m)
%          P.corrs_OAFlux  (pipelines/asi_regimes_parameter_calculation.m)
% Outputs: P.monte_carlo   synthetic correlations, threshold, and the
%                          fraction of observed curves below it

clear; clc

P = asi_regimes_paths;

nIterations = 1000;  % number of synthetic series
seed = 0;            % random seed, for reproducibility

%% Observed anomalies (same processing as the ASI regime calculation)

load(P.oaflux_NIO, 'flux_turb_NIO')

dt = datetime(datevec(flux_turb_NIO.time));
[~, keep_noleap] = doy_noleap(dt);

Qturb = flux_turb_NIO.lhf(:,:,keep_noleap) + flux_turb_NIO.shf(:,:,keep_noleap);
SST = flux_turb_NIO.sst(:,:,keep_noleap);

SST(SST > 50) = nan;
Qturb(isnan(SST)) = nan;

[SSTanom, SSTtendanom, Qanom, doy365] = asi_anomalies_noleap(SST, Qturb, dt(keep_noleap));
clear flux_turb_NIO SST Qturb

%% Randomized time series

nT = numel(doy365);

igood = ~isnan(Qanom) & ~isnan(SSTanom) & ~isnan(SSTtendanom);
Q_anom_good = Qanom(igood);
SST_anom_good = SSTanom(igood);
SST_tend_anom_good = SSTtendanom(igood);
clear SSTanom SSTtendanom Qanom

rng(seed, 'twister')
[Q_MC, SST_MC, SST_tend_MC] = deal(nan(nIterations,1,nT));
for nn = 1:nIterations
    Q_MC(nn,1,:) = datasample(Q_anom_good, nT, 'Replace', false);
    SST_MC(nn,1,:) = datasample(SST_anom_good, nT, 'Replace', false);
    SST_tend_MC(nn,1,:) = datasample(SST_tend_anom_good, nT, 'Replace', false);
end

%% Lagged correlations of the synthetic series (iteration x 1 x lag x day)

MC = lagged_corrs_doy(SST_MC, SST_tend_MC, Q_MC, doy365);

%% Noise threshold and fraction of observed curves below it (Eq. A1)

corr_mag_MC = squeeze(max(abs(MC.corrs_SST),[],3) + max(abs(MC.corrs_tend),[],3));
corr_thresh = prctile(corr_mag_MC(:), 95);

load(P.corrs_OAFlux, 'sst_thf_anomaly_corrs')
corr_mag_obs = squeeze(max(abs(sst_thf_anomaly_corrs.corrs_SST),[],3) + ...
                       max(abs(sst_thf_anomaly_corrs.corrs_tend),[],3));
frac_below = sum(corr_mag_obs(:) < corr_thresh) / sum(~isnan(corr_mag_obs(:)));

fprintf('95th percentile threshold of max|C_Ts,Q| + max|C_dTs/dt,Q|: %.3f\n', corr_thresh)
fprintf('Observed curves below threshold: %.2f%%\n', 100*frac_below)

%% Save

MC.nIterations = nIterations;
MC.seed = seed;
MC.corr_thresh_95 = corr_thresh;
MC.frac_obs_below_thresh = frac_below;

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.monte_carlo, 'MC', '-v7.3')
