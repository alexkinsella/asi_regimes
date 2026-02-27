% Make fig 3 of pre/during/post monsoon regimes and spatial variability
% Alex Kinsella, February 2026

load /Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/oaflux_NIO.mat
load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_OAFlux_260210.mat')

%% Calculations for maps
land_mask = double(isnan(asi_regimes.param(:,:,1)));
land_mask(land_mask == 0) = nan;

params_masked = asi_regimes.param;
params_masked_land = params_masked;
for dd = 1:size(params_masked,3)
    params_now = params_masked(:,:,dd);
    params_now(land_mask==1) = -999;
    params_masked_land(:,:,dd) = params_now;
end

%% Make curves for regions

% Build 2D lon/lat matrices
[lon2d, lat2d] = meshgrid(asi_regimes.lon, asi_regimes.lat);
lon2d = lon2d';   % match data orientation
lat2d = lat2d';

% Slanted Arabian Sea boundary
m = (62 - 50) / (26 - 5);
b = 50 - m * 5;

% 2D masks
reg1_mask2d = ...
    lon2d >= 50 & lon2d <= 80 & ...
    lat2d >= 5  & lat2d <= 26 & ...
    lon2d >= (m * lat2d + b);

reg2_mask2d = ...
    lon2d >= 80 & lon2d <= 98 & ...
    lat2d >= 5  & lat2d <= 23;

reg3_mask2d = ...
    lon2d >= 50 & lon2d <= 98 & ...
    lat2d >= -5 & lat2d <= 5;

% Compute regional averages
nt = size(params_masked,3);

reg1_meanparam = nan(nt,1);
reg2_meanparam = nan(nt,1);
reg3_meanparam = nan(nt,1);

for dd = 1:nt
    slice = params_masked(:,:,dd);

    reg1_meanparam(dd) = nanmean(slice(reg1_mask2d));
    reg2_meanparam(dd) = nanmean(slice(reg2_mask2d));
    reg3_meanparam(dd) = nanmean(slice(reg3_mask2d));
end
dates = datenum(2020,1,1):datenum(2020,1,365);
ndaysmean = 30;
reg1_movmean = periodic_movmean(reg1_meanparam,ndaysmean);
reg2_movmean = periodic_movmean(reg2_meanparam,ndaysmean);
reg3_movmean = periodic_movmean(reg3_meanparam,ndaysmean);

%% Plotting

dates = datetime(2020,1,1:365);
day1 = 136;
day2 = 183;
day3 = 306;

cm = cmocean('balance');
cm(1,:) = 0;

fig1 = makefig;
subplot(2,2,1)
pcolor(flux_turb_NIO.lon,flux_turb_NIO.lat,params_masked_land(:,:,day1)')
shading flat
clim([0 1]);
colormap(cm);
makespruce(30)
maketitle(['Pre-Monsoon: ',datestr(dates(day1),'mmm DD')],40)

subplot(2,2,2)
pcolor(flux_turb_NIO.lon,flux_turb_NIO.lat,params_masked_land(:,:,day2)')
shading flat
clim([0 1]);
colormap(cm);
makespruce(30)
maketitle(['During Monsoon: ',datestr(dates(day2),'mmm DD')],40)

subplot(2,2,3)
pcolor(flux_turb_NIO.lon,flux_turb_NIO.lat,params_masked_land(:,:,day3)')
shading flat
clim([0 1]);
colormap(cm);
makespruce(30)
maketitle(['Post-Monsoon: ',datestr(dates(day3),'mmm DD')],40)

subplot(2,2,4)
hold on
plot(dates,reg1_movmean,'linewidth',3)
plot(dates,reg2_movmean,'linewidth',3)
plot(dates,reg3_movmean,'linewidth',3)
yline(0.5,'linewidth',2,'color','k','linestyle','--')
hold off
leg = legend('Arabian Sea','Bay of Bengal','Equatorial','interpreter','latex');
xtickformat('MMM')
ax = gca;
ax.XTickLabel = ax.XTickLabel;
ylabel('ASI Parameter $P_o$','interpreter','latex')
makespruce(30)
set(leg,'fontsize',20)
maketitle('Regional ASI Regime Parameters',40)

% Create inset map
ax_inset = axes('Position', [0.6 0.3 0.1 0.1]);

set(ax_inset, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 2); % Add black border
hold on;
load coastlines;
plot(coastlon, coastlat, 'color', 'k', 'linewidth', 2);
reg3_lonmin = 50; reg3_lonmax = 98; reg3_latmin = -5; reg3_latmax = 5;
rectangle('Position', [reg3_lonmin, reg3_latmin, reg3_lonmax - reg3_lonmin, reg3_latmax - reg3_latmin], 'linewidth', 2,'edgecolor',"#EDB120");
% Arabian Sea polygon
AS_lon = [50 80 80 62 50];
AS_lat = [ 5  5 26 26  5];
% Bay of Bengal polygon
BoB_lon = [80 98 98 80 80];
BoB_lat = [ 5  5 23 23  5];
patch(AS_lon,  AS_lat,  'k','Edgecolor',"#0072BD", 'FaceColor','none', 'LineWidth',2);
patch(BoB_lon, BoB_lat, 'k','Edgecolor',"#D95319", 'FaceColor','none', 'LineWidth',2);
xlim([40 100]);
ylim([-5 30]);

sgtitle('Climatological Air-Sea Interaction Regimes','fontsize',50,'interpreter','latex')

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