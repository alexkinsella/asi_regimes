function [sst_thf_anomaly_corrs, asi_regimes] = asi_regimes_calculation(SST, Qturb, dt, lon, lat)
%ASI_REGIMES_CALCULATION Daily air-sea interaction regime parameter P_o.
%
% [sst_thf_anomaly_corrs, asi_regimes] = asi_regimes_calculation(SST, Qturb, dt, lon, lat)
%
% Implements Section 2b and the appendix of Kinsella and Mahadevan (2026):
%   1) Remove Feb. 29, form daily anomalies of SST, SST tendency, and
%      Qturb from an 11-day circular climatology smoothed with a 31-day
%      circular moving mean.
%   2) For each day of year, compute lagged correlations (lags -30 to 30
%      days) of SST and SST tendency with Qturb in a 31-day circular
%      window, concatenated across years.
%   3) Compute effective p values (Bretherton et al. 1999) and
%      Benjamini-Hochberg q values, and convert them to sigmoid weights.
%   4) Decompose the correlation curves into ocean- and
%      atmosphere-controlled parts by weighted nonnegative least squares.
%
% Inputs:
%   SST    lon x lat x time, deg C
%   Qturb  lon x lat x time, latent + sensible heat flux, W/m^2, positive
%          upward (ocean to atmosphere)
%   dt     datetime vector of daily times (Feb. 29 is removed here if
%          present; removing it beforehand saves memory)
%   lon, lat  grid vectors
%
% Outputs:
%   sst_thf_anomaly_corrs  lagged correlations and naive p values
%                          (lon x lat x lag x day of year)
%   asi_regimes            P_o as asi_regimes.param (lon x lat x day of year)
%
% Requires the Statistics and Machine Learning Toolbox (corr, tcdf).
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Alex Kinsella, June 2026

%% 1: 365-day no-leap calendar, tendency, climatology, and anomalies

[SSTanom, SSTtendanom, Qanom, doy365, dt_noleap, clim_opts] = asi_anomalies_noleap(SST, Qturb, dt);

nlon = numel(lon);
nlat = numel(lat);

%% 2: Lagged correlations for sliding day-of-year windows

sst_thf_anomaly_corrs = lagged_corrs_doy(SSTanom, SSTtendanom, Qanom, doy365);
sst_thf_anomaly_corrs.calendar = '365-day no-leap; Feb. 29 omitted; circular day-of-year windows';
sst_thf_anomaly_corrs.clim_halfwidth_days = clim_opts.clim_halfwidth_days;
sst_thf_anomaly_corrs.clim_smooth_days = clim_opts.clim_smooth_days;

corrs_SST = sst_thf_anomaly_corrs.corrs_SST;
corrs_tend = sst_thf_anomaly_corrs.corrs_tend;
lags = sst_thf_anomaly_corrs.lags;
nwind_half = sst_thf_anomaly_corrs.corr_halfwidth_days;

%% 3: Effective p values and Benjamini-Hochberg false discovery rate

% Nominal number of samples per day-of-year correlation window
nyears = numel(unique(year(dt_noleap)));
N = nyears * (2*nwind_half + 1);

% Lag-1 autocorrelations by grid point
[rho_SST,rho_tend,rho_Q] = deal(nan(nlon,nlat));

for ii = 1:nlon
    for jj = 1:nlat

        x1 = squeeze(SSTanom(ii,jj,1:end-1));
        x2 = squeeze(SSTanom(ii,jj,2:end));
        good = ~isnan(x1) & ~isnan(x2);
        if sum(good) > 3
            rho_SST(ii,jj) = corr(x1(good),x2(good));
        end

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

clear SSTanom SSTtendanom Qanom

N_eff_SST = N*(1-rho_SST.*rho_Q)./(1+rho_SST.*rho_Q);
N_eff_tend = N*(1-rho_tend.*rho_Q)./(1+rho_tend.*rho_Q);

q_SST_eff = effective_q_values(corrs_SST, N_eff_SST);
q_tend_eff = effective_q_values(corrs_tend, N_eff_tend);

%% 4: Sigmoid weights

k = 5; % Steepness parameter for sigmoid

weight_SST = 1-(1./(1+exp(-k*(2*q_SST_eff-1))));
weight_tend = 1-(1./(1+exp(-k*(2*q_tend_eff-1))));

clear q_SST_eff q_tend_eff

%% 5: Weighted nonnegative least squares decomposition

c = nan(nlon,nlat,365,2);

for dd = 1:365
    for ii = 1:nlon
        for jj = 1:nlat

            cSST = squeeze(corrs_SST(ii,jj,:,dd));
            cTend = squeeze(corrs_tend(ii,jj,:,dd));

            if all(isnan(cSST)) || all(isnan(cTend))
                continue
            end

            % Symmetric (even) and antisymmetric (odd) parts about zero lag
            SST_even = 1/2*(cSST+flip(cSST));
            SST_odd = 1/2*(cSST-flip(cSST));
            tend_even = 1/2*(cTend+flip(cTend));
            tend_odd = 1/2*(cTend-flip(cTend));

            if norm(SST_even)==0 || norm(SST_odd)==0 || ...
               norm(tend_even)==0 || norm(tend_odd)==0 || ...
               norm(cSST)==0 || norm(cTend)==0
                continue
            end

            SST_even_norm = SST_even/norm(SST_even);
            SST_odd_norm = SST_odd/norm(SST_odd);
            tend_even_norm = tend_even/norm(tend_even);
            tend_odd_norm = tend_odd/norm(tend_odd);

            SST_to_tend = norm(cSST)/norm(cTend);

            weights = diag([sqrt(squeeze(weight_SST(ii,jj,:,dd))); ...
                            sqrt(squeeze(weight_tend(ii,jj,:,dd)))]);

            % Column 1: ocean control (symmetric SST, antisymmetric tendency)
            % Column 2: atmosphere control (antisymmetric SST, symmetric tendency)
            C = weights*[[SST_even_norm; tend_odd_norm./SST_to_tend], ...
                         [SST_odd_norm;  tend_even_norm./SST_to_tend]];
            d = weights*[cSST;cTend];

            if any(isnan(C(:))) || any(isnan(d(:)))
                continue
            end

            c(ii,jj,dd,:) = lsqnonneg(C,d);

        end
    end
    fprintf('ASI regimes: done with day %d of 365\n', dd)
end

for ii = 1:nlon
    for jj = 1:nlat
        if norm(squeeze(c(ii,jj,1,:)))==0
            c(ii,jj,:,:) = nan;
        end
    end
end

asi_regimes.lon = lon;
asi_regimes.lat = lat;
asi_regimes.param = c(:,:,:,1).^2./(c(:,:,:,1).^2+c(:,:,:,2).^2);
asi_regimes.doy = 1:365;
asi_regimes.calendar = '365-day no-leap; Feb. 29 omitted; circular day-of-year windows';
asi_regimes.lags = lags;
asi_regimes.clim_halfwidth_days = clim_opts.clim_halfwidth_days;
asi_regimes.clim_smooth_days = clim_opts.clim_smooth_days;
asi_regimes.corr_halfwidth_days = nwind_half;

end

function q_eff = effective_q_values(corrs, N_eff)
%EFFECTIVE_Q_VALUES Benjamini-Hochberg q values from effective p values.

dof_eff = N_eff-2;
t_val_eff = corrs.*sqrt(dof_eff./(1-corrs.^2));
p_eff = 2*(1-tcdf(abs(t_val_eff),repmat(dof_eff,[1 1 size(corrs,[3 4])])));

idxGood = ~isnan(p_eff);
q_eff = nan(size(p_eff));
q_eff(idxGood) = bh_qvalues(p_eff(idxGood));
end
