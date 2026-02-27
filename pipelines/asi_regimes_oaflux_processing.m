% Process OAFlux files

lonmin = 40;
lonmax = 100;
latmin = -10;
latmax = 30;
region = 'NIO';

localpath = '/Volumes/proj/mahadevanlab/datasets/OAFlux/daily/';
lhf_files = dir([localpath,'turbulence/lh*.nc']);
shf_files = dir([localpath,'turbulence/sh*.nc']);
sst_files = dir([localpath,'turbulence/ts*.nc']);
qnet_files = dir([localpath,'netheat/qnet*.nc']);
lw_files = dir([localpath,'radiation/lw*.nc']);
sw_files = dir([localpath,'radiation/sw*.nc']);
qa_files = dir([localpath,'turbulence/qa*.nc']);
ta_files = dir([localpath,'turbulence/ta*.nc']);
ws_files = dir([localpath,'turbulence/ws*.nc']);
evap_files = dir([localpath,'evaporation/evapr*.nc']);

%% Process turb and evap files
lon = ncread([lhf_files(1).folder,'/',lhf_files(1).name],'lon');
lat = ncread([lhf_files(1).folder,'/',lhf_files(1).name],'lat');
loncheck = lon >= lonmin & lon <= lonmax;
latcheck = lat >= latmin & lat <= latmax;
reglon = lon(loncheck);
reglat = lat(latcheck);
time = [];
time_lens = [];
for yy = 1:numel(lhf_files)
    newtime = datenum(1985-1+yy,1,1) - 1 + ncread([lhf_files(yy).folder,'/',lhf_files(yy).name],'time');
    time = [time; newtime];
    time_lens = [time_lens; numel(newtime)];
    disp(['Done time year ',num2str(yy/numel(lhf_files))])
end
inds = [0; cumsum(time_lens)];

[lhf,shf,sst,qa,ta,ws,evap] = deal(nan(numel(reglon),numel(reglat),numel(time)));

for yy = 1:numel(lhf_files)
    try
        lhftmp = ncread([lhf_files(yy).folder,'/',lhf_files(yy).name],'lhtfl');
    catch 
        lhftmp = ncread([lhf_files(yy).folder,'/',lhf_files(yy).name],'lhflx');
    end
    try
        shftmp = ncread([shf_files(yy).folder,'/',shf_files(yy).name],'shtfl');
    catch
        shftmp = ncread([shf_files(yy).folder,'/',shf_files(yy).name],'shflx');
    end
    ssttmp = ncread([sst_files(yy).folder,'/',sst_files(yy).name],'tmpsf');
    qatmp = ncread([qa_files(yy).folder,'/',qa_files(yy).name],'hum2m');
    tatmp = ncread([ta_files(yy).folder,'/',ta_files(yy).name],'tmp2m');
    wstmp = ncread([ws_files(yy).folder,'/',ws_files(yy).name],'wnd10');
    evtmp = ncread([evap_files(yy).folder,'/',evap_files(yy).name],'evapr');

    lhf(:,:,inds(yy)+1:inds(yy+1)) = lhftmp(loncheck,latcheck,:);
    shf(:,:,inds(yy)+1:inds(yy+1)) = shftmp(loncheck,latcheck,:);
    sst(:,:,inds(yy)+1:inds(yy+1)) = ssttmp(loncheck,latcheck,:);
    qa(:,:,inds(yy)+1:inds(yy+1)) = qatmp(loncheck,latcheck,:);
    ta(:,:,inds(yy)+1:inds(yy+1)) = tatmp(loncheck,latcheck,:);
    ws(:,:,inds(yy)+1:inds(yy+1)) = wstmp(loncheck,latcheck,:);
    evap(:,:,inds(yy)+1:inds(yy+1)) = evtmp(loncheck,latcheck,:);

    disp(['Done with turbulence year ',num2str(yy),' of ',num2str(numel(lhf_files))])
end
lhf(lhf>3e3) = nan;
shf(shf>3e3) = nan;
sst(sst>3e2) = nan;
qa(qa>3e2) = nan;
ta(ta>3e2) = nan;
ws(ws>3e2) = nan;
evap(evap>3e3) = nan;

%% Save
time = double(time);
flux_turb_NIO.lhf = lhf; 
flux_turb_NIO.shf = shf; 
flux_turb_NIO.sst = sst;
flux_turb_NIO.ta = ta;
flux_turb_NIO.qa = qa;
flux_turb_NIO.ws = ws;
flux_turb_NIO.evap = evap;
flux_turb_NIO.lon = reglon;
flux_turb_NIO.lat = reglat;
flux_turb_NIO.time = time;
%save('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/oaflux_NIO.mat','flux_turb_NIO')
