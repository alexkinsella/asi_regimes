% Process ERA5 turbulent heat fluxes and SST over the Indian Ocean into
% daily means, 1985-2020.
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  ERA5 6-hourly NetCDF files, one per year
%          (P.era5_dir/ERA5_THF_SST_IO_<year>_rates.nc) with variables
%          msshf and mslhf (mean surface sensible and latent heat flux,
%          W/m^2, ECMWF convention: positive downward) and sst (K)
% Outputs: P.era5_1985_2020, structure ERA5_THF_SST_IO with fields
%          lon, lat, dnum (datenum), dvec, shf, lhf (W/m^2, positive
%          downward), and sst (deg C)

clear; clc

P = asi_regimes_paths;

yrs = 1985:2020;

%% Grid

f1 = fullfile(P.era5_dir, sprintf('ERA5_THF_SST_IO_%d_rates.nc', yrs(1)));
lat = ncread(f1,'latitude');
lon = ncread(f1,'longitude');

dnum_all = datenum(yrs(1),1,1):datenum(yrs(end),12,31);
[shf_daily,lhf_daily,sst_daily] = deal(nan(numel(lon),numel(lat),numel(dnum_all)));

%% Daily means of 6-hourly data, one year at a time

cnt = 1;
for yy = 1:numel(yrs)
    filename = fullfile(P.era5_dir, sprintf('ERA5_THF_SST_IO_%d_rates.nc', yrs(yy)));

    time = ncread(filename,'time');
    dnum = datenum(datetime(1900,1,1) + hours(time));
    ndays = numel(dnum)/4;

    shf = ncread(filename,'msshf');
    lhf = ncread(filename,'mslhf');
    sst = ncread(filename,'sst');

    for dy = 1:ndays
        curdate = datenum(yrs(yy),1,1)+dy-1;
        inds = dnum>=curdate & dnum<curdate+1;
        shf_daily(:,:,cnt+dy-1) = mean(shf(:,:,inds),3,'omitnan');
        lhf_daily(:,:,cnt+dy-1) = mean(lhf(:,:,inds),3,'omitnan');
        sst_daily(:,:,cnt+dy-1) = mean(sst(:,:,inds),3,'omitnan');
    end
    cnt = cnt+ndays;

    fprintf('Done with year %d\n', yrs(yy))
end

%% Save

ERA5_THF_SST_IO.lon = lon;
ERA5_THF_SST_IO.lat = lat;
ERA5_THF_SST_IO.dnum = dnum_all;
ERA5_THF_SST_IO.dvec = datevec(dnum_all);
ERA5_THF_SST_IO.shf = shf_daily;
ERA5_THF_SST_IO.lhf = lhf_daily;
ERA5_THF_SST_IO.sst = sst_daily-273.15;

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.era5_1985_2020,'ERA5_THF_SST_IO','-v7.3')
