% Influence of radiation on the ASI regimes (Section 4b): distribution of
% simultaneous correlations between SST tendency anomalies and net surface
% radiative flux anomalies, which estimate the proportion of model
% stochastic forcing contributed by radiation (Eq. 16).
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Anomalies and lagged correlations use the same 365-day no-leap,
% circular day-of-year method as the ASI regime calculation, with CERES
% net radiation (SW + LW, positive downward) in place of Qturb, over the
% CERES-OAFlux overlap (March 2000 - December 2020). Correlations are
% kept where the Benjamini-Hochberg q value of the effective p value is
% at most 0.1.
%
% Inputs:  P.ceres_rad  (pipelines/asi_regimes_ceres_processing.m)
%          P.oaflux_NIO (pipelines/asi_regimes_oaflux_processing.m), SST
% Outputs: printed summary of the correlation distribution;
%          P.fig_dir/asi_regimes_radiative_correlation.png

clear; clc

P = asi_regimes_paths;

alpha = 0.1;  % false discovery rate

%% Load and align SST and radiative flux on the common grid and dates

load(P.ceres_rad, 'CERES_rad')
load(P.oaflux_NIO, 'flux_turb_NIO')

assert(isequal(double(CERES_rad.lon(:)), double(flux_turb_NIO.lon(:))) && ...
       isequal(double(CERES_rad.lat(:)), double(flux_turb_NIO.lat(:))), ...
       'CERES and OAFlux grids differ')

ceres_days = floor(CERES_rad.dnum(:));
[dnum, iC, iO] = intersect(ceres_days, flux_turb_NIO.time(:));

Qrad = CERES_rad.SWnet(:,:,iC) + CERES_rad.LWnet(:,:,iC);
SST = flux_turb_NIO.sst(:,:,iO);
clear CERES_rad flux_turb_NIO

SST(SST > 50) = nan;
Qrad(isnan(SST)) = nan;

%% Anomalies and lagged correlations (same method as the ASI regimes)

[SSTanom, SSTtendanom, Qanom, doy365, dt_noleap] = ...
    asi_anomalies_noleap(SST, Qrad, datetime(dnum,'ConvertFrom','datenum'));
clear SST Qrad

C = lagged_corrs_doy(SSTanom, SSTtendanom, Qanom, doy365);

%% Effective p values and Benjamini-Hochberg q values for the
%% SST tendency-radiation correlations

nyears = numel(unique(year(dt_noleap)));
N = nyears * (2*C.corr_halfwidth_days + 1);

% Lag-1 autocorrelations by grid point
[rho_tend, rho_Q] = deal(nan(size(Qanom,[1 2])));
for ii = 1:size(Qanom,1)
    for jj = 1:size(Qanom,2)
        x1 = squeeze(SSTtendanom(ii,jj,1:end-1));
        x2 = squeeze(SSTtendanom(ii,jj,2:end));
        good = ~isnan(x1) & ~isnan(x2);
        if sum(good) > 3
            rho_tend(ii,jj) = corr(x1(good),x2(good));
        end

        x1 = squeeze(Qanom(ii,jj,1:end-1));
        x2 = squeeze(Qanom(ii,jj,2:end));
        good = ~isnan(x1) & ~isnan(x2);
        if sum(good) > 3
            rho_Q(ii,jj) = corr(x1(good),x2(good));
        end
    end
end

N_eff_tend = N*(1-rho_tend.*rho_Q)./(1+rho_tend.*rho_Q);
dof_eff_tend = N_eff_tend-2;

corrs_tend = C.corrs_tend;
t_val_tend_eff = corrs_tend.*sqrt(dof_eff_tend./(1-corrs_tend.^2));
p_tend_eff = 2*(1-tcdf(abs(t_val_tend_eff),repmat(dof_eff_tend,[1 1 size(corrs_tend,[3 4])])));

idxGood = ~isnan(p_tend_eff);
q_tend_eff = nan(size(p_tend_eff));
q_tend_eff(idxGood) = bh_qvalues(p_tend_eff(idxGood));

%% Simultaneous (zero-lag) correlations that pass the FDR test

ilag0 = find(C.lags == 0);
tend_rad_sim_corr = squeeze(corrs_tend(:,:,ilag0,:));
q_tend_eff_sim = squeeze(q_tend_eff(:,:,ilag0,:));
tend_rad_sim_corr(q_tend_eff_sim > alpha) = nan;

r = tend_rad_sim_corr(~isnan(tend_rad_sim_corr));
fprintf('Significant simultaneous dSST/dt-Qrad correlations: %d of %d (%.1f%%)\n', ...
    numel(r), sum(~isnan(q_tend_eff_sim(:))), 100*numel(r)/sum(~isnan(q_tend_eff_sim(:))))
fprintf('Median %.3f; fraction above 0.2: %.1f%%; above 0.3: %.1f%%\n', ...
    median(r), 100*mean(r > 0.2), 100*mean(r > 0.3))

%% Histogram

fig1 = makefig(1800, 987);
h = histogram(r);
[~, imax] = max(h.Values);
fprintf('Histogram peak at %.3f to %.3f\n', h.BinEdges(imax), h.BinEdges(imax+1))
xlabel('Correlation','interpreter','latex')
ylabel('Count','interpreter','latex')
hold on
xline(0,'--','color','k','linewidth',2)
makespruce(30)
maketitle('Simultaneous SST Tendency - Radiative Heat Flux Anomaly Correlation',40)
xlim([-0.3 0.4])

export_figure(fig1, P.fig_dir, 'asi_regimes_radiative_correlation')
