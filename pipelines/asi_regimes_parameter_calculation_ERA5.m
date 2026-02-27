% Calculate ASI regime parameters from ERA5

%% Load and process
load('/Volumes/proj/mahadevanlab/kinsella/data/asi_regimes/ERA5_THF_SST_IO_1985_2020.mat');

dvec = ERA5_THF_SST_IO.dvec;
Qturb = -ERA5_THF_SST_IO.lhf + -ERA5_THF_SST_IO.shf;
SST = ERA5_THF_SST_IO.sst;
nlon = numel(ERA5_THF_SST_IO.lon);
nlat = numel(ERA5_THF_SST_IO.lat);
SSTtend = nan(size(SST));
for ii = 1:size(SSTtend,1)
    for jj = 1:size(SSTtend,2)
        SSTtend(ii,jj,:) = gradient(squeeze(SST(ii,jj,:))); % K/day
    end
end

%% Create daily climo and anom
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
for tt = 1:numel(ERA5_THF_SST_IO.dnum)
    dd = day(datetime(dvec(tt,:)),'dayofyear');
    SSTanom(:,:,tt) = SST(:,:,tt) - SSTclim(:,:,dd);
    SSTtendanom(:,:,tt) = SSTtend(:,:,tt) - SSTtendclim(:,:,dd);
    Qanom(:,:,tt) = Qturb(:,:,tt) - Qclim(:,:,dd);
end

%% 1: Calculate lagged corr for sliding months
nlag = 30;
corrs_SST = nan(nlon,nlat,2*nlag+1,365);
pvals_SST = nan(nlon,nlat,2*nlag+1,365);
corrs_tend = nan(nlon,nlat,2*nlag+1,365);
pvals_tend = nan(nlon,nlat,2*nlag+1,365);

doys = day(datetime(dvec),'dayofyear');
nwind = 30; % Days in sliding window

for Qday = 1:365 % Choose fixed central day for Qturb
    corrs_SST_day = nan(nlon,nlat,2*nlag+1);
    corrs_tend_day = nan(nlon,nlat,2*nlag+1);
    pvals_SST_day = nan(nlon,nlat,2*nlag+1);
    pvals_tend_day = nan(nlon,nlat,2*nlag+1);
    Qfixed = Qanom(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
    for lag = -nlag:nlag % For positive lag, SST lags Qturb
        if lag >= 0
            SSTshift = circshift(SSTanom,-lag,3);
            tendshift = circshift(SSTtendanom,-lag,3);
            SSTseries = SSTshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            tendseries = tendshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            Qfixed_trim = Qfixed(:,:,1:end-lag);
            SSTseries = SSTseries(:,:,1:end-lag);
            tendseries = tendseries(:,:,1:end-lag);
        elseif lag < 0
            SSTshift = circshift(SSTanom,-lag,3);
            tendshift = circshift(SSTtendanom,-lag,3);
            SSTseries = SSTshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            tendseries = tendshift(:,:,doys>=Qday-nwind/2 & doys<=Qday+nwind/2);
            Qfixed_trim = Qfixed(:,:,-lag+1:end); 
            SSTseries = SSTseries(:,:,-lag+1:end);
            tendseries = tendseries(:,:,-lag+1:end);
        end
        for lon = 1:nlon
            for lat = 1:nlat
                [A,P] = corrcoef(squeeze(SSTseries(lon,lat,:)),squeeze(Qfixed_trim(lon,lat,:)));
                corrs_SST_day(lon,lat,lag+nlag+1) = A(2,1);
                pvals_SST_day(lon,lat,lag+nlag+1) = P(2,1);
                
                [A,P] = corrcoef(squeeze(tendseries(lon,lat,:)),squeeze(Qfixed_trim(lon,lat,:)));
                corrs_tend_day(lon,lat,lag+nlag+1) = A(2,1);
                pvals_tend_day(lon,lat,lag+nlag+1) = P(2,1);
            end
        end
        disp(['Done with lag ',num2str(lag)])
    end
    corrs_SST(:,:,:,Qday) = corrs_SST_day;
    pvals_SST(:,:,:,Qday) = pvals_SST_day;
    corrs_tend(:,:,:,Qday) = corrs_tend_day;
    pvals_tend(:,:,:,Qday) = pvals_tend_day;
    disp(['Done with day ',num2str(Qday)])
end

sst_thf_anomaly_corrs.corrs_SST = corrs_SST;
sst_thf_anomaly_corrs.corrs_tend = corrs_tend;
sst_thf_anomaly_corrs.pvals_SST = pvals_SST;
sst_thf_anomaly_corrs.pvals_tend = pvals_tend;
%save('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/sst_thf_anomaly_corrs_ERA5_260211.mat','sst_thf_anomaly_corrs','-v7.3')

%% Load if resuming here
load('/Volumes/proj/mahadevanlab/kinsella/data/asi_regimes/ERA5_THF_SST_IO_1985_2020.mat');

dvec = ERA5_THF_SST_IO.dvec;
Qturb = -ERA5_THF_SST_IO.lhf + -ERA5_THF_SST_IO.shf;
SST = ERA5_THF_SST_IO.sst;
SSTtend = nan(size(SST));
for ii = 1:size(SSTtend,1)
    for jj = 1:size(SSTtend,2)
        SSTtend(ii,jj,:) = gradient(squeeze(SST(ii,jj,:))); % K/day
    end
end

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
for tt = 1:numel(ERA5_THF_SST_IO.dnum)
    dd = day(datetime(dvec(tt,:)),'dayofyear');
    SSTanom(:,:,tt) = SST(:,:,tt) - SSTclim(:,:,dd);
    SSTtendanom(:,:,tt) = SSTtend(:,:,tt) - SSTtendclim(:,:,dd);
    Qanom(:,:,tt) = Qturb(:,:,tt) - Qclim(:,:,dd);
end

clear ERA5_THF_SST_IO dvec nhalf nsmooth SST Qturb SSTtend Qclim SSTclim SSTtendclim 

load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/sst_thf_anomaly_corrs_ERA5_260211.mat')
corrs_SST = sst_thf_anomaly_corrs.corrs_SST;
corrs_tend = sst_thf_anomaly_corrs.corrs_tend;
%lon = sst_thf_anomaly_corrs.lon;
%lat = sst_thf_anomaly_corrs.lat;

clear sst_thf_anomaly_corrs

% 2: Calculate p-eff and do FDR test 
% Number of time points 
N = 36*30; % 30 days per year for 36 years
% Find lag-1 autocorrelations by grid point
[rho_SST,rho_tend,rho_Q] = deal(nan(size(corrs_SST,[1 2])));
for ii = 1:size(corrs_SST,1)
    for jj = 1:size(corrs_SST,2)
        rho_SST(ii,jj) = corr(squeeze(SSTanom(ii,jj,1:end-1)),squeeze(SSTanom(ii,jj,2:end)));
        rho_tend(ii,jj) = corr(squeeze(SSTtendanom(ii,jj,1:end-1)),squeeze(SSTtendanom(ii,jj,2:end)));
        rho_Q(ii,jj) = corr(squeeze(Qanom(ii,jj,1:end-1)),squeeze(Qanom(ii,jj,2:end)));
    end
    disp(['Done with ',num2str(ii),' of ',num2str(size(corrs_SST,1))])
end
N_eff_SST = N*(1-rho_SST.*rho_Q)./(1+rho_SST.*rho_Q);
N_eff_tend = N*(1-rho_tend.*rho_Q)./(1+rho_tend.*rho_Q);
clear rho_SST rho_tend rho_Q Qanom SSTanom SSTtendanom

% Calculate t value and effective p values for SST correlations
%dof = N-2;
dof_eff_SST = N_eff_SST-2;
%t_val_SST_naive = corrs_SST.*sqrt(dof./(1-corrs_SST.^2));
t_val_SST_eff = corrs_SST.*sqrt(dof_eff_SST./(1-corrs_SST.^2));
%p_SST_naive = 2*(1-tcdf(abs(t_val_SST_naive),dof));
p_SST_eff = 2*(1-tcdf(abs(t_val_SST_eff),repmat(dof_eff_SST,[1 1 61 365])));
clear t_val_SST_eff dof_eff_SST N_eff_SST

% Apply Benjamini-Hochberg for SST correlations
idxGood = ~isnan(p_SST_eff);
pvals_good = p_SST_eff(idxGood);
qvals_good = mafdr(pvals_good, 'BHFDR', true);
clear pvals_good
q_SST_eff = nan(size(p_SST_eff));
clear p_SST_eff
q_SST_eff(idxGood) = qvals_good;
clear qvals_good idxGood

% Calculate t value and effective p values for tendency correlations
%dof = N-2;
dof_eff_tend = N_eff_tend-2;
%t_val_tend_naive = corrs_tend.*sqrt(dof./(1-corrs_tend.^2));
t_val_tend_eff = corrs_tend.*sqrt(dof_eff_tend./(1-corrs_tend.^2));
%p_tend_naive = 2*(1-tcdf(abs(t_val_tend_naive),dof));
p_tend_eff = 2*(1-tcdf(abs(t_val_tend_eff),repmat(dof_eff_tend,[1 1 61 365])));
clear t_val_tend_eff dof_eff_tend N_eff_tend

% Apply Benjamini-Hochberg for tendency correlations
idxGood = ~isnan(p_tend_eff);
pvals_good = p_tend_eff(idxGood);
qvals_good = mafdr(pvals_good, 'BHFDR', true);
clear pvals_good
q_tend_eff = nan(size(p_tend_eff));
clear p_tend_eff
q_tend_eff(idxGood) = qvals_good;
clear qvals_good idxGood

% 3: Create sigmoid weighting functions
k = 5; % Steepness parameter for sigmoid
%alpha = 0.05; % FDR parameter
weight_SST = 1-(1./(1+exp(-k*(2*q_SST_eff-1))));
weight_tend = 1-(1./(1+exp(-k*(2*q_tend_eff-1))));
clear q_SST_eff
clear q_tend_eff

%% 4: Calculate ASI regimes
lags = -30:30;
c = nan([size(squeeze(corrs_tend(:,:,1,:))),2]);
resnorm = nan(size(squeeze(corrs_tend(:,:,1,:))));
residual = nan([size(squeeze(corrs_tend(:,:,1,:))),2*numel(lags)]);

for dd = 1:size(corrs_tend,4)
    for ii = 1:size(corrs_tend,1)
        for jj = 1:size(corrs_tend,2)

            SST = squeeze(corrs_SST(ii,jj,:,dd));
            tend = squeeze(corrs_tend(ii,jj,:,dd));

            SST_even = 1/2*(SST+flip(SST));
            SST_odd = 1/2*(SST-flip(SST));
            tend_even = 1/2*(tend+flip(tend));
            tend_odd = 1/2*(tend-flip(tend));

            SST_even_norm = SST_even/norm(SST_even);
            SST_odd_norm = SST_odd/norm(SST_odd);
            tend_even_norm = tend_even/norm(tend_even);
            tend_odd_norm = tend_odd/norm(tend_odd);

            SST_to_tend = norm(SST)/norm(tend);

            weights = diag([sqrt(squeeze(weight_SST(ii,jj,:,dd)));sqrt(squeeze(weight_tend(ii,jj,:,dd)))]);

            C = weights*[[SST_even_norm; tend_odd_norm./SST_to_tend], [SST_odd_norm; tend_even_norm./SST_to_tend]];
            d = weights*[SST;tend];

            [c(ii,jj,dd,:),resnorm(ii,jj,dd),residual(ii,jj,dd,:)] = lsqnonneg(C,d);
        
        end
    end
    disp(['Done with day ',num2str(dd)])
end

for ii = 1:size(c,1)
    for jj = 1:size(c,2)
        if norm(squeeze(c(ii,jj,1,:)))==0
            c(ii,jj,:,:) = nan;
        end
    end
end

asi_regimes.lon = ERA5_THF_SST_IO.lon;
asi_regimes.lat = ERA5_THF_SST_IO.lat;
asi_regimes.param = c(:,:,:,1).^2./(c(:,:,:,1).^2+c(:,:,:,2).^2);
save('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_ERA5_260220.mat','asi_regimes','-v7.3')
