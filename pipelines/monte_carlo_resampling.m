%% Test significance of lagged correlations based on a Monte Carlo resampling method
% Alex Kinsella, February 2026

load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/sst_thf_anomaly_corrs_260209.mat')
load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/flux_turb_NIO_ext.mat')

dvec = datevec(flux_turb_NIO.time);
Qturb = flux_turb_NIO.lhf + flux_turb_NIO.shf;
SST = flux_turb_NIO.sst;
SST(SST>50) = nan;
Qturb(isnan(SST)) = nan;
nlon = numel(flux_turb_NIO.lon);
nlat = numel(flux_turb_NIO.lat);
SSTtend = nan(size(SST));
for ii = 1:size(SSTtend,1)
    for jj = 1:size(SSTtend,2)
        SSTtend(ii,jj,:) = gradient(squeeze(SST(ii,jj,:))); % K/day
    end
end

% Create daily climo and anom
nhalf = 5;
nsmooth = 30;
SSTclim = clim_daily(SST,dvec,nhalf);
SSTtendclim = clim_daily(SSTtend,dvec,nhalf);
Qclim = clim_daily(Qturb,dvec,nhalf);
SSTclim = movmean(SSTclim,nsmooth,3);
SSTtendclim = movmean(SSTtendclim,nsmooth,3);
Qclim = movmean(Qclim,nsmooth,3);

SSTanom = nan(size(SST));
SSTtendanom = nan(size(SSTtend));
Qanom = nan(size(Qturb));
for tt = 1:numel(flux_turb_NIO.time)
    dd = day(datetime(dvec(tt,:)),'dayofyear');
    SSTanom(:,:,tt) = SST(:,:,tt) - SSTclim(:,:,dd);
    SSTtendanom(:,:,tt) = SSTtend(:,:,tt) - SSTtendclim(:,:,dd);
    Qanom(:,:,tt) = Qturb(:,:,tt) - Qclim(:,:,dd);
end

% Create randomized time series
nIterations = 1000; % Number of iterations for MC
nT = size(SSTanom,3); % Length of time series

doys = repmat(1:365,[1 ceil(nT/365)]);
doys = doys(1:nT);

igood = ~isnan(Qanom);
Q_anom_good = Qanom(igood);
SST_anom_good = SSTanom(igood);
SST_tend_anom_good = SSTtendanom(igood);

[Q_MC_samples,SST_MC_samples,...
    SST_tend_MC_samples] = deal(nan(1,nIterations,nT));
for nn = 1:nIterations
    Q_MC_samples(1,nn,:) = datasample(Q_anom_good,nT,'Replace',false);
    SST_MC_samples(1,nn,:) = datasample(SST_anom_good,nT,'Replace',false);
    SST_tend_MC_samples(1,nn,:) = datasample(SST_tend_anom_good,nT,'Replace',false);
end

% Constants
nlag = 30; % Number of lags to compute 
nwind = 30; % Days in sliding window

corrs_SST_MC = nan(1, nIterations, 2*nlag+1, 365);
pvals_SST_MC = nan(1, nIterations, 2*nlag+1, 365);
corrs_tend_MC = nan(1, nIterations, 2*nlag+1, 365);
pvals_tend_MC = nan(1, nIterations, 2*nlag+1, 365);

for Qday = 1:365 % Choose fixed central day for Qturb
    [corrs_SST_day,pvals_SST_day,corrs_tend_day,pvals_tend_day] = deal(nan(1,nIterations,2*nlag+1));
    Qfixed = Q_MC_samples(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
    for lag = -nlag:nlag % For positive lag, SST lags Qturb
        if lag >= 0
            SSTshift = circshift(SST_MC_samples,-lag,3);
            tendshift = circshift(SST_tend_MC_samples,-lag,3);
            SSTseries = SSTshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            tendseries = tendshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            Qfixed_trim = Qfixed(:,:,1:end-lag);
            SSTseries = SSTseries(:,:,1:end-lag);
            tendseries = tendseries(:,:,1:end-lag);
        elseif lag < 0
            SSTshift = circshift(SST_MC_samples,-lag,3);
            tendshift = circshift(SST_tend_MC_samples,-lag,3);
            SSTseries = SSTshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            tendseries = tendshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            Qfixed_trim = Qfixed(:,:,-lag+1:end); 
            SSTseries = SSTseries(:,:,-lag+1:end);
            tendseries = tendseries(:,:,-lag+1:end);
        end
        for nn = 1:nIterations
                [A,P] = corrcoef(squeeze(SSTseries(1,nn,:)),squeeze(Qfixed_trim(1,nn,:)));
                corrs_SST_day(1,nn,lag+nlag+1) = A(2,1);
                pvals_SST_day(1,nn,lag+nlag+1) = P(2,1);
                
                [A,P] = corrcoef(squeeze(tendseries(1,nn,:)),squeeze(Qfixed_trim(1,nn,:)));
                corrs_tend_day(1,nn,lag+nlag+1) = A(2,1);
                pvals_tend_day(1,nn,lag+nlag+1) = P(2,1);
            
        end
        disp(['Done with lag ',num2str(lag)])
    end
    corrs_SST_MC(:,:,:,Qday) = corrs_SST_day;
    pvals_SST_MC(:,:,:,Qday) = pvals_SST_day;
    corrs_tend_MC(:,:,:,Qday) = corrs_tend_day;
    pvals_tend_MC(:,:,:,Qday) = pvals_tend_day;
    disp(['Done with day ',num2str(Qday)])
end

MC.corrs_SST = corrs_SST_MC;
MC.corrs_tend = corrs_tend_MC;
MC.pvals_SST = pvals_SST_MC;
MC.pvals_tend = pvals_tend_MC;

%save('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/sst_thf_anomaly_corrs_monte_carlo_260210.mat','MC','-v7.3')

% Find distribution of correlation amplitudes
corr_mag_MC = squeeze(max(abs(corrs_SST_MC),[],3)+max(abs(corrs_tend_MC),[],3));
% find 95 percentile
corr_thresh_bootstrap = prctile(corr_mag_MC(:),95);

% Create mask for obs
corr_mag_obs = squeeze(max(abs(sst_thf_anomaly_corrs.corrs_SST),[],3)+max(abs(sst_thf_anomaly_corrs.corrs_tend),[],3));
sum(corr_mag_obs(:)<0.194)/sum(~isnan(corr_mag_obs(:)))
