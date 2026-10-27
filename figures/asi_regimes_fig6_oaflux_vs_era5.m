% Figure 6: ERA5 vs. OAFlux ASI regimes
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Also makes supplemental Video 1 (daily OAFlux and ERA5 ASI regimes) when
% make_video = true, and prints the overall OAFlux-ERA5 correlation
% (Section 3c).
%
% Inputs:  P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
%          P.asi_ERA5   (pipelines/asi_regimes_parameter_calculation_ERA5.m)
% Outputs: P.fig_dir/asi_regimes_fig6_oaflux_vs_era5.png
%          P.fig_dir/asi_regimes_video1_oaflux_vs_era5.mp4 (if make_video)

clear; clc

P = asi_regimes_paths;

make_video = false;  % set true to write supplemental Video 1 (365 frames)

% Load ASI regimes
asi_era5 = load(P.asi_ERA5, 'asi_regimes');
asi_OA = load(P.asi_OAFlux, 'asi_regimes');

lonmatrix_era5 = repmat(asi_era5.asi_regimes.lon,[1 size(asi_era5.asi_regimes.param,[2 3])]);
latmatrix_era5 = permute(repmat(asi_era5.asi_regimes.lat,[1 size(asi_era5.asi_regimes.param,[1 3])]),[2 1 3]);

lonmatrix_OA = repmat(asi_OA.asi_regimes.lon,[1 size(asi_OA.asi_regimes.param,[2 3])]);
latmatrix_OA = permute(repmat(asi_OA.asi_regimes.lat,[1 size(asi_OA.asi_regimes.param,[1 3])]),[2 1 3]);

% Calculations for maps

oa_params_masked_land = land_mask(asi_OA.asi_regimes.param);
era5_params_masked_land = land_mask(asi_era5.asi_regimes.param);

% Difference maps: ERA5 - OAFlux, interpolated onto OAFlux grid
% Positive difference = ERA5 has larger P_o, i.e. more ocean control.
lon_era5 = asi_era5.asi_regimes.lon(:);
lat_era5 = asi_era5.asi_regimes.lat(:);

% griddedInterpolant requires monotonic coordinate vectors
[lon_era5_sort, ilon] = sort(lon_era5);
[lat_era5_sort, ilat] = sort(lat_era5);

era5_param_sort = asi_era5.asi_regimes.param(ilon, ilat, :);

era5_on_OA = nan(size(asi_OA.asi_regimes.param));

F = griddedInterpolant( ...
    {lon_era5_sort, lat_era5_sort}, ...
    era5_param_sort(:,:,1), ...
    'linear', ...
    'none');

for dd = 1:size(asi_OA.asi_regimes.param,3)
    F.Values = era5_param_sort(:,:,dd);
    era5_on_OA(:,:,dd) = F(lonmatrix_OA(:,:,1), latmatrix_OA(:,:,1));
end

diff_era5_minus_OA = era5_on_OA - asi_OA.asi_regimes.param;

cm = cmocean('balance');
cm(1,:) = 0;

cm_diff = cmocean('diff');
diff_lim = 1;

% -------------------------------------------------------------------------
% Subset plotted fields to visible map domain
% This avoids rendering global pcolor objects and then cropping visually.
% -------------------------------------------------------------------------
plot_lonlim = [40 100];
plot_latlim = [-5 30];

lon_OA_vec   = asi_OA.asi_regimes.lon(:);
lat_OA_vec   = asi_OA.asi_regimes.lat(:);
lon_era5_vec = asi_era5.asi_regimes.lon(:);
lat_era5_vec = asi_era5.asi_regimes.lat(:);

iOA = lon_OA_vec >= plot_lonlim(1) & lon_OA_vec <= plot_lonlim(2);
jOA = lat_OA_vec >= plot_latlim(1) & lat_OA_vec <= plot_latlim(2);

iERA5 = lon_era5_vec >= plot_lonlim(1) & lon_era5_vec <= plot_lonlim(2);
jERA5 = lat_era5_vec >= plot_latlim(1) & lat_era5_vec <= plot_latlim(2);

lonplot_OA = lonmatrix_OA(iOA,jOA,1);
latplot_OA = latmatrix_OA(iOA,jOA,1);

lonplot_era5 = lonmatrix_era5(iERA5,jERA5,1);
latplot_era5 = latmatrix_era5(iERA5,jERA5,1);

oa_plot   = oa_params_masked_land(iOA,jOA,:);
era5_plot = era5_params_masked_land(iERA5,jERA5,:);
diff_plot = diff_era5_minus_OA(iOA,jOA,:);

%% Make side by side plot for manuscript figure

load('coastlines.mat');
coast_lon = coastlon;
coast_lat = coastlat;

fig1 = makefig(1800, 987);  % window size of the published figure
set(fig1, 'Renderer', 'opengl')

% -------------------------------------------------------------------------
% May 15
% -------------------------------------------------------------------------
doy = 135; % May 15

ax = subplot(2,3,1);
pcolor(lonplot_OA, latplot_OA, oa_plot(:,:,doy));
shading flat
colormap(ax, cm)
clim([0 1])
makespruce(20)
hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', 0.75)
hold off
maketitle('(a) May 15: OAFlux',30)
xlim([40 100])
ylim([-5 30])
ylabel('Latitude ($^\circ$N)', 'Interpreter','latex','fontsize',20)

ax = subplot(2,3,2);
pcolor(lonplot_era5, latplot_era5, era5_plot(:,:,doy));
shading flat
colormap(ax, cm)
clim([0 1])
makespruce(20)
hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', 0.75)
hold off
maketitle('(b) May 15: ERA5',30)
xlim([40 100])
ylim([-5 30])

ax = subplot(2,3,3);
pcolor(lonplot_OA, latplot_OA, diff_plot(:,:,doy));
shading flat
colormap(ax, cm_diff)
clim([-diff_lim diff_lim])
set(ax, 'Color', 'k')  % show NaNs / land as black
makespruce(20)
hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', 0.75)
hold off
maketitle('(c) May 15: ERA5 - OAFlux',30)
xlim([40 100])
ylim([-5 30])

% -------------------------------------------------------------------------
% July 1
% -------------------------------------------------------------------------
doy = 182; % July 1

ax = subplot(2,3,4);
pcolor(lonplot_OA, latplot_OA, oa_plot(:,:,doy));
shading flat
colormap(ax, cm)
clim([0 1])
makespruce(20)
hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', 0.75)
hold off
maketitle('(d) July 1: OAFlux',30)
xlim([40 100])
ylim([-5 30])
ylabel('Latitude ($^\circ$N)', 'Interpreter','latex','fontsize',20)
xlabel('Longitude ($^\circ$E)', 'Interpreter','latex','fontsize',20)

ax = subplot(2,3,5);
pcolor(lonplot_era5, latplot_era5, era5_plot(:,:,doy));
shading flat
colormap(ax, cm)
clim([0 1])
makespruce(20)
hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', 0.75)
hold off
maketitle('(e) July 1: ERA5',30)
xlim([40 100])
ylim([-5 30])
xlabel('Longitude ($^\circ$E)', 'Interpreter','latex','fontsize',20)

ax = subplot(2,3,6);
pcolor(lonplot_OA, latplot_OA, diff_plot(:,:,doy));
shading flat
colormap(ax, cm_diff)
clim([-diff_lim diff_lim])
set(ax, 'Color', 'k')  % show NaNs / land as black
makespruce(20)
hold on
plot(coast_lon, coast_lat, 'w', 'LineWidth', 0.75)
hold off
maketitle('(f) July 1: ERA5 - OAFlux',30)
xlim([40 100])
ylim([-5 30])
xlabel('Longitude ($^\circ$E)', 'Interpreter','latex','fontsize',20)

% -------------------------------------------------------------------------
% Main ASI colorbar on left
% -------------------------------------------------------------------------
ax_cb1 = axes('visible','off');
colormap(ax_cb1, cm)
clim(ax_cb1, [0 1])

cb = makecb([0 1], cm, 'ASI Parameter $P_o$', 20, ax_cb1);
pos = get(cb,'position');
set(cb,'position',pos + [-.82 -0.02 0 0])
cb.Label.FontSize = 25;

% y=1 label
text( ...
    -0.16, 0.93, sprintf('%s\n%s','Ocean','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',20, ...
    'Interpreter','latex','rotation', 90)

% y=0 label
text( ...
    -0.12, -0.05, sprintf('%s\n%s','Atm.','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',20, ...
    'Interpreter','latex','Rotation', 90)

% -------------------------------------------------------------------------
% Difference colorbar on right
% -------------------------------------------------------------------------
ax_cb2 = axes('visible','off');
colormap(ax_cb2, cm_diff)
clim(ax_cb2, [-diff_lim diff_lim])

cb_diff = makecb( ...
    [-diff_lim diff_lim], ...
    cm_diff, ...
    '$\Delta P_o$ (ERA5 - OAFlux)', ...
    20, ...
    ax_cb2);
cb_diff.Label.FontSize = 25;

pos = get(cb_diff,'position');
set(cb_diff,'position',pos + [0.05 -0.02 0 0])

sgtitle('ASI Regimes in OAFlux vs. ERA5','interpreter','latex','fontsize',40)

export_figure(fig1, P.fig_dir, 'asi_regimes_fig6_oaflux_vs_era5')

%% Supplemental Video 1: daily ASI regimes, side by side

if make_video
    if ~exist(P.fig_dir,'dir')
        mkdir(P.fig_dir)
    end
    if any(strcmp({VideoWriter.getProfiles.Name}, 'MPEG-4'))
        vw = VideoWriter(fullfile(P.fig_dir, 'asi_regimes_video1_oaflux_vs_era5'), 'MPEG-4');
    else  % MPEG-4 is not available on Linux
        vw = VideoWriter(fullfile(P.fig_dir, 'asi_regimes_video1_oaflux_vs_era5'), 'Motion JPEG AVI');
    end
    vw.FrameRate = 10;
    open(vw)

    % Day of year in the 365-day no-leap calendar; 2023 is a non-leap year
    dates = datetime(2023,1,1:365);

    fig2 = makefig(1800, 600);
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

        % Redrawn every frame: subplot deletes the overlapping colorbar axes
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

        sgtitle(['ASI Regimes on ',datestr(dates(dd),'mmm dd')],'interpreter','latex','fontsize',50)

        writeVideo(vw, print(fig2, '-RGBImage', '-r100'))
    end
    close(vw)
end

%% Overall OAFlux-ERA5 correlation (all grid points and all days)

OA_corr   = asi_OA.asi_regimes.param;
ERA5_corr = era5_on_OA;

% Restrict to the common valid ocean points
good = isfinite(OA_corr) & isfinite(ERA5_corr);

OA_vec   = OA_corr(good);
ERA5_vec = ERA5_corr(good);

R = corr(OA_vec, ERA5_vec, 'Type','Pearson');

RMSE = sqrt(mean((ERA5_vec - OA_vec).^2));
bias = mean(ERA5_vec - OA_vec);

fprintf('\nOverall OAFlux-ERA5 agreement\n')
fprintf('Number of gridpoint-days: %d\n', numel(OA_vec))
fprintf('Pearson R = %.4f\n', R)
fprintf('R^2 = %.4f\n', R^2)
fprintf('RMSE = %.4f\n', RMSE)
fprintf('Mean bias (ERA5 - OAFlux) = %.4f\n\n', bias)

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