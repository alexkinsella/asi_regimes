% ERA5 vs. OAFlux ASI regimes
% Alex Kinsella, February 2026

% Load ASI regimes
asi_era5= load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_ERA5_260220.mat');
asi_OA = load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_OAFlux_260210.mat');

% Subset to NIO region
% Bounding box
lonmin = 50;
lonmax = 100;
latmin = 5;
latmax = 26;

% Line from (5N, 50E) to (26N, 62E)
m = (62 - 50) / (26 - 5);  % slope
b = 50 - m * 5;            % intercept

lonmatrix_era5 = repmat(asi_era5.asi_regimes.lon,[1 size(asi_era5.asi_regimes.param,[2 3])]);
latmatrix_era5 = permute(repmat(asi_era5.asi_regimes.lat,[1 size(asi_era5.asi_regimes.param,[1 3])]),[2 1 3]);

lonmatrix_OA = repmat(asi_OA.asi_regimes.lon,[1 size(asi_OA.asi_regimes.param,[2 3])]);
latmatrix_OA = permute(repmat(asi_OA.asi_regimes.lat,[1 size(asi_OA.asi_regimes.param,[1 3])]),[2 1 3]);

era5_check = (lonmatrix_era5 >= lonmin & lonmatrix_era5 <= lonmax & ...
latmatrix_era5 >= latmin & latmatrix_era5 <= latmax & ...
lonmatrix_era5 >= (m * latmatrix_era5 + b));

OA_check = (lonmatrix_OA >= lonmin & lonmatrix_OA <= lonmax & ...
latmatrix_OA >= latmin & latmatrix_OA <= latmax & ...
lonmatrix_OA >= (m * latmatrix_OA + b));

param_masked_era5 = asi_era5.asi_regimes.param;
param_masked_era5(~era5_check) = nan;

param_masked_OA = asi_OA.asi_regimes.param;
param_masked_OA(~OA_check) = nan;

%% Make side by side plot for a date

% Calculations for maps

oa_params_masked_land = land_mask(asi_OA.asi_regimes.param);
era5_params_masked_land = land_mask(asi_era5.asi_regimes.param);

cm = cmocean('balance');
cm(1,:) = 0;

load coastlines


fig1 = makefig; 

doy = 135; % May 15

subplot(2,2,1)
pcolor(lonmatrix_OA(:,:,1),latmatrix_OA(:,:,1),oa_params_masked_land(:,:,doy));
shading flat
%cb = makecb([0 1],cmocean('balance'),'ASI Regime',30);
colormap(cm)
clim([0 1])
makespruce(30)
maketitle('May 15: OAFlux',40)
xlim([40 100])
ylim([-5 30])
%daspect([1 cos(pi/180*15) 1])

subplot(2,2,2)
pcolor(lonmatrix_era5(:,:,1),latmatrix_era5(:,:,1),era5_params_masked_land(:,:,doy));
shading flat
colormap(cm)
clim([0 1])
makespruce(30)
maketitle('May 15: ERA5',40)
xlim([40 100])
ylim([-5 30])
%daspect([1 cos(pi/180*15) 1])

doy = 182; % July 1

subplot(2,2,3)
pcolor(lonmatrix_OA(:,:,1),latmatrix_OA(:,:,1),oa_params_masked_land(:,:,doy));
shading flat
%cb = makecb([0 1],cmocean('balance'),'ASI Regime',30);
colormap(cm)
clim([0 1])
makespruce(30)
maketitle('July 1: OAFlux',40)
xlim([40 100])
ylim([-5 30])
%daspect([1 cos(pi/180*15) 1])

subplot(2,2,4)
pcolor(lonmatrix_era5(:,:,1),latmatrix_era5(:,:,1),era5_params_masked_land(:,:,doy));
shading flat
colormap(cm)
clim([0 1])
makespruce(30)
maketitle('July 1: ERA5',40)
xlim([40 100])
ylim([-5 30])
%daspect([1 cos(pi/180*15) 1])

ax2 = axes('visible','off');
cb = makecb([0 1],cmocean('balance'),'ASI Parameter $P_o$',30,ax2);
pos = get(cb,'position');
set(cb,'position',pos + [-.83 -0.02 0 0])

% y=1 label
text( ...
    -0.16, 0.93, sprintf('%s\n%s','Ocean','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',25, ...
    'Interpreter','latex','rotation', 90)

% y=0 label
text( ...
    -0.12, -0.05, sprintf('%s\n%s','Atm.','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',25, ...
    'Interpreter','latex','Rotation', 90)

sgtitle('ASI Regimes in OAFlux vs. ERA5','interpreter','latex','fontsize',50)

%% Make side-by-side video

load coastlines

giffing = 1;
fig1 = makefig(1800, 600); 
for dd = 1:size(asi_OA.asi_regimes.param,3)
    subplot(1,2,1)
    pcolor(lonmatrix_OA(:,:,1),latmatrix_OA(:,:,1),oa_params_masked_land(:,:,dd));
    shading flat
    colormap(cm)
    clim([0 1])
    makespruce(30)
    maketitle('OAFlux',40)
    xlim([40 100])
    ylim([-5 30])
    daspect([1 cos(pi/180*15) 1])
    
    subplot(1,2,2)
    pcolor(lonmatrix_era5(:,:,1),latmatrix_era5(:,:,1),era5_params_masked_land(:,:,dd));
    shading flat
    colormap(cm)
    clim([0 1])
    makespruce(30)
    maketitle('ERA5',40)
    xlim([40 100])
    ylim([-5 30])
    daspect([1 cos(pi/180*15) 1])
    
    ax2 = axes('visible','off');
    cb = makecb([0 1],cmocean('balance'),'ASI Parameter $P_o$',30,ax2);
    pos = get(cb,'position');
    set(cb,'position',pos + [-.83 -0.02 0 0])
    
    % y=1 label
    text( ...
        -0.16, 0.93, sprintf('%s\n%s','Ocean','control'), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize',25, ...
        'Interpreter','latex','rotation', 90)
    
    % y=0 label
    text( ...
        -0.12, -0.05, sprintf('%s\n%s','Atm.','control'), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom', ...
        'FontSize',25, ...
        'Interpreter','latex','Rotation', 90)

    sgtitle(['ASI Regimes on ',datestr(datetime(2023,1,dd),'mmm dd')],'interpreter','latex','fontsize',50)

    if giffing==1 && dd==1
        gif('/Users/akinsella/Projects/asi_regimes/figures/final/asi_regimes_video1_OA_vs_ERA5.gif','DelayTime',0.1)
    elseif giffing==1 && dd~=1
        gif
    elseif giffing==0
        pause(0.01)
    end
    pause(0.1)
end

function params_masked_land = land_mask(params)

land_mask = double(isnan(params(:,:,1)));
land_mask(land_mask == 0) = nan;

params_masked = params;
params_masked_land = params_masked;
for dd = 1:size(params_masked,3)
    params_now = params_masked(:,:,dd);
    params_now(land_mask==1) = -999;
    params_masked_land(:,:,dd) = params_now;
end

end