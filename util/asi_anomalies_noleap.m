function [SSTanom, SSTtendanom, Qanom, doy365, dt_noleap, clim_opts] = asi_anomalies_noleap(SST, Qturb, dt)
%ASI_ANOMALIES_NOLEAP Daily anomalies of SST, SST tendency, and Qturb.
%
% [SSTanom, SSTtendanom, Qanom, doy365, dt_noleap, clim_opts] = asi_anomalies_noleap(SST, Qturb, dt)
%
% Removes Feb. 29, computes the SST tendency as a centered difference,
% and subtracts a daily climatology: an 11-day circular day-of-year mean
% (Hobday et al. 2016) smoothed with a 31-day circular moving mean.
%
% Inputs:
%   SST, Qturb  lon x lat x time
%   dt          datetime vector of daily times
%
% Outputs:
%   SSTanom, SSTtendanom, Qanom  lon x lat x time, Feb. 29 removed
%   doy365      365-day no-leap day of year for each time
%   dt_noleap   times with Feb. 29 removed
%   clim_opts   climatology window settings
%
% Alex Kinsella, June 2026

[doy365, keep_noleap, dt_noleap] = doy_noleap(dt);
if ~all(keep_noleap)
    SST = SST(:,:,keep_noleap);
    Qturb = Qturb(:,:,keep_noleap);
end

nlon = size(SST,1);
nlat = size(SST,2);
nt = numel(doy365);

% SST tendency on no-leap daily time axis
SSTtend = nan(size(SST));
for ii = 1:nlon
    for jj = 1:nlat
        SSTtend(ii,jj,:) = gradient(squeeze(SST(ii,jj,:))); % K/day
    end
end

nhalf = 5;          % 11-day climatology window: +/- 5 days
nsmooth = 31;       % 31-day circular smoothing window

SSTclim = clim_daily_noleap(SST, doy365, nhalf);
SSTtendclim = clim_daily_noleap(SSTtend, doy365, nhalf);
Qclim = clim_daily_noleap(Qturb, doy365, nhalf);

SSTclim = circular_movmean_doy(SSTclim, nsmooth);
SSTtendclim = circular_movmean_doy(SSTtendclim, nsmooth);
Qclim = circular_movmean_doy(Qclim, nsmooth);

SSTanom = nan(size(SST));
SSTtendanom = nan(size(SSTtend));
Qanom = nan(size(Qturb));

for tt = 1:nt
    dd = doy365(tt);
    SSTanom(:,:,tt) = SST(:,:,tt) - SSTclim(:,:,dd);
    SSTtendanom(:,:,tt) = SSTtend(:,:,tt) - SSTtendclim(:,:,dd);
    Qanom(:,:,tt) = Qturb(:,:,tt) - Qclim(:,:,dd);
end

clim_opts.clim_halfwidth_days = nhalf;
clim_opts.clim_smooth_days = nsmooth;

end
