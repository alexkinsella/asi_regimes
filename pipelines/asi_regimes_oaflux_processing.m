% Process daily OAFlux turbulent fluxes, SST, and near-surface variables
% over the northern Indian Ocean, 1985-2020.
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  OAFlux daily NetCDF files, one per year and variable
%          (P.oaflux_dir/turbulence/ and P.oaflux_dir/evaporation/)
% Outputs: P.oaflux_NIO, structure flux_turb_NIO with fields
%          lhf, shf (W/m^2), sst, ta (deg C), qa (g/kg), ws (m/s),
%          evap (cm/yr), lon, lat, and time (datenum)

clear; clc

P = asi_regimes_paths;

%% Region and period

lonmin = 40;
lonmax = 100;
latmin = -10;
latmax = 30;

yrs = 1985:2020;

%% Grid and time

f1 = fullfile(P.oaflux_dir, 'turbulence', sprintf('lh_oaflux_%d.nc', yrs(1)));
lon = ncread(f1,'lon');
lat = ncread(f1,'lat');
loncheck = lon >= lonmin & lon <= lonmax;
latcheck = lat >= latmin & lat <= latmax;
reglon = lon(loncheck);
reglat = lat(latcheck);

time = [];
time_lens = nan(numel(yrs),1);
for yy = 1:numel(yrs)
    f = fullfile(P.oaflux_dir, 'turbulence', sprintf('lh_oaflux_%d.nc', yrs(yy)));
    newtime = datenum(yrs(yy),1,1) - 1 + ncread(f,'time');
    time = [time; newtime]; %#ok<AGROW>
    time_lens(yy) = numel(newtime);
end
inds = [0; cumsum(time_lens)];

%% Read and subset each year

[lhf,shf,sst,qa,ta,ws,evap] = deal(nan(numel(reglon),numel(reglat),numel(time)));

for yy = 1:numel(yrs)
    fname = @(subdir,prefix) fullfile(P.oaflux_dir, subdir, sprintf('%s_oaflux_%d.nc', prefix, yrs(yy)));

    % Latent and sensible heat flux variable names differ between years
    try
        lhftmp = ncread(fname('turbulence','lh'),'lhtfl');
    catch
        lhftmp = ncread(fname('turbulence','lh'),'lhflx');
    end
    try
        shftmp = ncread(fname('turbulence','sh'),'shtfl');
    catch
        shftmp = ncread(fname('turbulence','sh'),'shflx');
    end
    ssttmp = ncread(fname('turbulence','ts'),'tmpsf');
    qatmp  = ncread(fname('turbulence','qa'),'hum2m');
    tatmp  = ncread(fname('turbulence','ta'),'tmp2m');
    wstmp  = ncread(fname('turbulence','ws'),'wnd10');
    evtmp  = ncread(fname('evaporation','evapr'),'evapr');

    tt = inds(yy)+1:inds(yy+1);
    lhf(:,:,tt)  = lhftmp(loncheck,latcheck,:);
    shf(:,:,tt)  = shftmp(loncheck,latcheck,:);
    sst(:,:,tt)  = ssttmp(loncheck,latcheck,:);
    qa(:,:,tt)   = qatmp(loncheck,latcheck,:);
    ta(:,:,tt)   = tatmp(loncheck,latcheck,:);
    ws(:,:,tt)   = wstmp(loncheck,latcheck,:);
    evap(:,:,tt) = evtmp(loncheck,latcheck,:);

    fprintf('Done with year %d\n', yrs(yy))
end

% Remove fill values
lhf(lhf>3e3) = nan;
shf(shf>3e3) = nan;
sst(sst>3e2) = nan;
qa(qa>3e2) = nan;
ta(ta>3e2) = nan;
ws(ws>3e2) = nan;
evap(evap>3e3) = nan;

%% Save

flux_turb_NIO.lhf = lhf;
flux_turb_NIO.shf = shf;
flux_turb_NIO.sst = sst;
flux_turb_NIO.ta = ta;
flux_turb_NIO.qa = qa;
flux_turb_NIO.ws = ws;
flux_turb_NIO.evap = evap;
flux_turb_NIO.lon = reglon;
flux_turb_NIO.lat = reglat;
flux_turb_NIO.time = double(time);

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.oaflux_NIO,'flux_turb_NIO')
