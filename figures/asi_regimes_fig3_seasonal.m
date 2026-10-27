% Figure 3: Pre/during/post monsoon regimes and spatial variability
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs: P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
% Output: P.fig_dir/asi_regimes_fig3_seasonal.png

clear; clc

P = asi_regimes_paths;

load(P.asi_OAFlux, 'asi_regimes')

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

    reg1_meanparam(dd) = mean(slice(reg1_mask2d),'omitnan');
    reg2_meanparam(dd) = mean(slice(reg2_mask2d),'omitnan');
    reg3_meanparam(dd) = mean(slice(reg3_mask2d),'omitnan');
end
ndaysmean = 30;
reg1_movmean = periodic_movmean(reg1_meanparam,ndaysmean);
reg2_movmean = periodic_movmean(reg2_meanparam,ndaysmean);
reg3_movmean = periodic_movmean(reg3_meanparam,ndaysmean);

%% Plotting

S = shaperead(P.coast_10m, 'UseGeoCoords', true);

coast_lon = [S.Lon];
coast_lat = [S.Lat];

% Day of year in the 365-day no-leap calendar; 2021 is a non-leap year,
% so dates(dd) gives the matching calendar date
dates = datetime(2021,1,1:365);
day1 = 136;   % May 16
day2 = 183;   % Jul 2
day3 = 306;   % Nov 2

cm = cmocean('balance');
cm(1,:) = 0;

fig1 = makefig(1800, 987);  % window size of the published figure
clf

% -------------------------------------------------------------------------
% Manual layout
% -------------------------------------------------------------------------
% Normalized figure positions: [left bottom width height]

% Use more vertical space: taller panels, top row higher, bottom row lower
pos_map1 = [0.145 0.535 0.350 0.325];
pos_map2 = [0.565 0.535 0.350 0.325];
pos_map3 = [0.145 0.095 0.350 0.325];

% Time-series panel: match panel c height/bottom
pos_ts   = [0.565 0.095 0.295 0.325];

map_xlim = [40 100];
map_ylim = [-5 30];

map_font   = 20;
title_font = 30;
coast_lw   = 0.55;

% Mooring locations
moor_lon = [61.5 89.5];
moor_lat = [15.5 18.0];
moor_marker_size = 20;

% -------------------------------------------------------------------------
% Panel a: pre-monsoon
% -------------------------------------------------------------------------
ax1 = axes('Position', pos_map1);
pcolor(asi_regimes.lon, asi_regimes.lat, params_masked_land(:,:,day1)')
shading flat
clim([0 1]);
colormap(cm);

hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', coast_lw)
plot(moor_lon, moor_lat, 'p', ...
    'MarkerSize', moor_marker_size, ...
    'MarkerFaceColor', 'y', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.6)
hold off

xlim(map_xlim)
ylim(map_ylim)
makespruce(map_font)
ylabel('Latitude ($^\circ$N)', 'Interpreter','latex', 'FontSize', map_font)
maketitle(['(a) Pre-Monsoon: ', datestr(dates(day1),'mmm DD')], title_font)

% -------------------------------------------------------------------------
% Panel b: during monsoon
% -------------------------------------------------------------------------
ax2 = axes('Position', pos_map2);
pcolor(asi_regimes.lon, asi_regimes.lat, params_masked_land(:,:,day2)')
shading flat
clim([0 1]);
colormap(cm);

hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', coast_lw)
plot(moor_lon, moor_lat, 'p', ...
    'MarkerSize', moor_marker_size, ...
    'MarkerFaceColor', 'y', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.6)
hold off

xlim(map_xlim)
ylim(map_ylim)
makespruce(map_font)
% No y-label here; left-column label is enough
maketitle(['(b) During Monsoon: ', datestr(dates(day2),'mmm DD')], title_font)

% -------------------------------------------------------------------------
% Panel c: post-monsoon
% -------------------------------------------------------------------------
ax3 = axes('Position', pos_map3);
pcolor(asi_regimes.lon, asi_regimes.lat, params_masked_land(:,:,day3)')
shading flat
clim([0 1]);
colormap(cm);

hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', coast_lw)
plot(moor_lon, moor_lat, 'p', ...
    'MarkerSize', moor_marker_size, ...
    'MarkerFaceColor', 'y', ...
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 0.6)
hold off

xlim(map_xlim)
ylim(map_ylim)
makespruce(map_font)
xlabel('Longitude ($^\circ$E)', 'Interpreter','latex', 'FontSize', map_font)
ylabel('Latitude ($^\circ$N)', 'Interpreter','latex', 'FontSize', map_font)
maketitle(['(c) Post-Monsoon: ', datestr(dates(day3),'mmm DD')], title_font)

% -------------------------------------------------------------------------
% Panel d: regional curves
% -------------------------------------------------------------------------
ax4 = axes('Position', pos_ts);
hold on
plot(dates, reg1_movmean, 'linewidth', 3)
plot(dates, reg2_movmean, 'linewidth', 3)
plot(dates, reg3_movmean, 'linewidth', 3)
yline(0.5, 'linewidth', 2, 'color', 'k', 'linestyle', '--')
hold off

ax4.XLim = [datetime(2021,1,1) datetime(2022,1,1)];
ax4.XTick = [datetime(2021,1,1), datetime(2021,4,1), ...
             datetime(2021,7,1), datetime(2021,10,1), ...
             datetime(2022,1,1)];
ax4.XTickLabel = {'Jan','Apr','Jul','Oct','Jan'};
ax4.XTickLabelRotation = 0;
xlabel('')
ylabel('ASI Parameter $P_o$', 'interpreter','latex')
makespruce(map_font)
maketitle('(d) Regional ASI Regime Parameters', title_font)

% Put legend in the extra right-side space, not on top of the data
leg = legend('Arabian Sea','Bay of Bengal','Equatorial', ...
    'interpreter','latex', ...
    'Box','on');
set(leg, 'FontSize', 16, 'Units','normalized')
set(leg, 'Position', [0.865 0.295 0.120 0.095])

% -------------------------------------------------------------------------
% Inset map in panel d
% -------------------------------------------------------------------------
ax_inset = axes('Position', [0.585 0.340 0.100 0.090]);

set(ax_inset, 'Box', 'on', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1.5);
hold on;

% Use the same high-resolution coastline here
plot(coast_lon, coast_lat, 'color', 'k', 'linewidth', 1.0);

reg3_lonmin = 50; reg3_lonmax = 98; reg3_latmin = -5; reg3_latmax = 5;
rectangle('Position', ...
    [reg3_lonmin, reg3_latmin, reg3_lonmax - reg3_lonmin, reg3_latmax - reg3_latmin], ...
    'linewidth', 1.5, 'edgecolor', "#EDB120");

% Arabian Sea polygon
AS_lon = [50 80 80 62 50];
AS_lat = [ 5  5 26 26  5];

% Bay of Bengal polygon
BoB_lon = [80 98 98 80 80];
BoB_lat = [ 5  5 23 23  5];

patch(AS_lon,  AS_lat,  'k', 'Edgecolor', "#0072BD", ...
    'FaceColor','none', 'LineWidth', 1.5);
patch(BoB_lon, BoB_lat, 'k', 'Edgecolor', "#D95319", ...
    'FaceColor','none', 'LineWidth', 1.5);

xlim([40 100]);
ylim([-5 30]);
set(gca, 'XTick', [], 'YTick', [])
hold off

% -------------------------------------------------------------------------
% Title
% -------------------------------------------------------------------------
sgtitle('Climatological Air-Sea Interaction Regimes', ...
    'fontsize', 44, 'interpreter','latex')

% -------------------------------------------------------------------------
% Colorbar
% -------------------------------------------------------------------------

% Remove any old colorbars if re-running this section interactively
delete(findall(gcf, 'Type', 'ColorBar'))

% Dummy invisible axes just to host the shared colorbar
ax_cb = axes('Position', [0.01 0.01 0.01 0.01], 'Visible', 'off');
colormap(ax_cb, cm)
clim(ax_cb, [0 1])

cb = colorbar(ax_cb);
set(cb, 'Position', [0.050 0.170 0.015 0.660])
set(cb, 'FontSize', 20)
set(cb, 'TickLabelInterpreter', 'latex')
cb.Ticks = 0:0.1:1;
cb.TickLabels = arrayfun(@(x) sprintf('$%g$', x), cb.Ticks, ...
    'UniformOutput', false);
cb.Label.String = '';

% Colorbar labels
ax_lab = axes('Position', [0 0 1 1], 'Visible', 'off');

% Main colorbar label, left of the colorbar
text(ax_lab, 0.027, 0.500, 'ASI Parameter $P_o$', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',25, ...
    'Interpreter','latex', ...
    'Rotation',90)

% y=1 label
text(ax_lab, 0.015, 0.800, sprintf('%s\n%s','Ocean','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',20, ...
    'Interpreter','latex', ...
    'Rotation',90)

% y=0 label
text(ax_lab, 0.015, 0.200, sprintf('%s\n%s','Atm.','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',20, ...
    'Interpreter','latex', ...
    'Rotation',90)

export_figure(fig1, P.fig_dir, 'asi_regimes_fig3_seasonal')
