function C = lagged_corrs_doy(SSTanom, SSTtendanom, Qanom, doy365)
%LAGGED_CORRS_DOY Lagged SST-Qturb and dSST/dt-Qturb correlations by day of year.
%
% C = lagged_corrs_doy(SSTanom, SSTtendanom, Qanom, doy365)
%
% For each day of year, Qturb anomalies are taken in a 31-day circular
% day-of-year window across all years, and correlated with SST and SST
% tendency anomalies shifted by lags of -30 to 30 days (positive lag:
% SST lags Qturb), Eq. (1) of Kinsella and Mahadevan (2026).
%
% Inputs:
%   SSTanom, SSTtendanom, Qanom  n1 x n2 x time (e.g. lon x lat x time)
%   doy365      365-day no-leap day of year for each time
%
% Output: structure C with fields
%   corrs_SST, corrs_tend, pvals_SST, pvals_tend   n1 x n2 x lag x day of year
%   lags, doy, corr_halfwidth_days
%
% Alex Kinsella, June 2026

n1 = size(Qanom,1);
n2 = size(Qanom,2);
nt = numel(doy365);

nlag = 30;
lags = -nlag:nlag;
nwind_half = 15;    % 31-day correlation window: +/- 15 days

corrs_SST = nan(n1,n2,2*nlag+1,365);
pvals_SST = nan(n1,n2,2*nlag+1,365);
corrs_tend = nan(n1,n2,2*nlag+1,365);
pvals_tend = nan(n1,n2,2*nlag+1,365);

for central_doy = 1:365

    corrs_SST_day = nan(n1,n2,2*nlag+1);
    corrs_tend_day = nan(n1,n2,2*nlag+1);
    pvals_SST_day = nan(n1,n2,2*nlag+1);
    pvals_tend_day = nan(n1,n2,2*nlag+1);

    % Times t for turbulent heat flux, selected by circular DOY window
    tind0 = find(circular_doy_window(doy365, central_doy, nwind_half));

    for ilag = 1:numel(lags)

        % Positive lag means SST/SST tendency occurs after Qturb
        tind = tind0;
        xind = tind + lags(ilag);

        valid = xind >= 1 & xind <= nt;
        tind = tind(valid);
        xind = xind(valid);

        Qfixed = Qanom(:,:,tind);
        SSTseries = SSTanom(:,:,xind);
        tendseries = SSTtendanom(:,:,xind);

        for ii = 1:n1
            for jj = 1:n2

                x = squeeze(SSTseries(ii,jj,:));
                y = squeeze(Qfixed(ii,jj,:));
                good = ~isnan(x) & ~isnan(y);

                if sum(good) > 3
                    [A,P] = corrcoef(x(good), y(good));
                    corrs_SST_day(ii,jj,ilag) = A(2,1);
                    pvals_SST_day(ii,jj,ilag) = P(2,1);
                end

                x = squeeze(tendseries(ii,jj,:));
                good = ~isnan(x) & ~isnan(y);

                if sum(good) > 3
                    [A,P] = corrcoef(x(good), y(good));
                    corrs_tend_day(ii,jj,ilag) = A(2,1);
                    pvals_tend_day(ii,jj,ilag) = P(2,1);
                end

            end
        end
    end

    corrs_SST(:,:,:,central_doy) = corrs_SST_day;
    pvals_SST(:,:,:,central_doy) = pvals_SST_day;
    corrs_tend(:,:,:,central_doy) = corrs_tend_day;
    pvals_tend(:,:,:,central_doy) = pvals_tend_day;

    fprintf('Correlations: done with day %d of 365\n', central_doy)
end

C.corrs_SST = corrs_SST;
C.corrs_tend = corrs_tend;
C.pvals_SST = pvals_SST;
C.pvals_tend = pvals_tend;
C.lags = lags;
C.doy = 1:365;
C.corr_halfwidth_days = nwind_half;

end
