% Figure 5: Relationship between ASI regimes and wind
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs: P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
%         P.oaflux_NIO (pipelines/asi_regimes_oaflux_processing.m)
% Output: P.fig_dir/asi_regimes_fig5_wind.png

clear; clc

P = asi_regimes_paths;

%% Load data
load(P.asi_OAFlux, 'asi_regimes')
load(P.oaflux_NIO, 'flux_turb_NIO')

%% Make OAFlux wind clim
% As in the published figure, the wind climatology uses the full
% (leap-year) calendar with a 30-day moving mean, and its first 365 days
% are paired with the 365-day no-leap P_o. After Feb. 28 the two are
% offset by one day; this has a negligible effect on the results.
nhalf = 5;
nsmooth = 30;

ws_clim_OAflux = clim_daily(flux_turb_NIO.ws, datevec(flux_turb_NIO.time), nhalf);
ws_clim_OAflux = movmean(ws_clim_OAflux, nsmooth, 3);
ws_clim_OAflux_short = ws_clim_OAflux(:,:,1:365);

%% Combine Arabian Sea and Bay of Bengal

lonmatrix = repmat(asi_regimes.lon,[1 size(asi_regimes.param,[2 3])]);
latmatrix = permute(repmat(asi_regimes.lat,[1 size(asi_regimes.param,[1 3])]),[2 1 3]);

% Line from (5N, 50E) to (26N, 62E)
m = (62 - 50) / (26 - 5);  % slope
b = 50 - m * 5;            % intercept

reg1_mask = ...
    lonmatrix >= 50 & lonmatrix <= 80 & ...
    latmatrix >= 5  & latmatrix <= 26 & ...
    lonmatrix >= (m * latmatrix + b);  % Arabian Sea mask

% Define spatial mask for Bay of Bengal
reg2_mask = ...
    lonmatrix >= 80 & lonmatrix <= 98 & ...
    latmatrix >= 5  & latmatrix <= 23;

% Define AS + BoB polygon outlines for plotting

% Arabian Sea polygon (rectangle clipped by the line lon = m*lat + b)
AS_lon = [50 80 80 62 50];
AS_lat = [ 5  5 26 26  5];

% Bay of Bengal rectangle
BoB_lon = [80 98 98 80 80];
BoB_lat = [ 5  5 23 23  5];

% Apply masks to 3D arrays
ws_both = ws_clim_OAflux_short(reg1_mask|reg2_mask);
asi_both = asi_regimes.param(reg1_mask|reg2_mask);

x_both = ws_both(:);
y_both = asi_both(:);

% Create 2D histogram
numBins = 100;  % Set the number of bins
[counts_both, edgesX_both, edgesY_both] = histcounts2(x_both, y_both, numBins);

% Define bin centers for plotting
binCentersX_both = (edgesX_both(1:end-1) + edgesX_both(2:end)) / 2;
binCentersY_both = (edgesY_both(1:end-1) + edgesY_both(2:end)) / 2;

% Median Y value in each X bin
medianYInXBin_both = nan(1, numBins);

for i = 1:numBins
    inXBin_both = x_both >= edgesX_both(i) & x_both < edgesX_both(i+1);
    if any(inXBin_both)
        medianYInXBin_both(i) = median(y_both(inXBin_both),'omitnan');
    end
end

%% Plotting
S = shaperead(P.coast_10m, 'UseGeoCoords', true);

% Plot the density as a 2D image
fig1=makefig(1800, 987);  % window size of the published figure

ax1 = subplot(2,2,3);
imagesc(binCentersX_both, binCentersY_both, (counts_both.'));
axis xy;
makecb([0 100],cmocean('amp'),'Count',20,ax1);
xlabel('Windspeed','interpreter','latex');
ylabel('ASI Parameter $P_o$','interpreter','latex');
xlim([2.7 12])
makespruce(20)
maketitle('(c) Windspeed-ASI Regime Density Plot',30);

% Overlay median Y values for each X bin as red dots
hold on;
for i = 1:numBins
    if ~isnan(medianYInXBin_both(i))
        plot(binCentersX_both(i), medianYInXBin_both(i), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    end
end
hold off;

% Now make maps

cm = cmocean('balance');
cm(1,:) = 0;

dd = 182;   % July 1 in the 365-day no-leap calendar

% Add panels
ax3 = subplot(2,2,2);
pcolor(asi_regimes.lon,asi_regimes.lat,asi_regimes.param(:,:,dd)')
shading flat
makecb([0 1],cm,'ASI Parameter $P_o$',20,ax3);
hold on
geoshow(S, ...
    'Color', 'w', ...
    'LineWidth', 1)
hold off
makespruce(20)
xlabel('Longitude ($^\circ$E)', 'Interpreter','latex','fontsize',20)
maketitle('(b) Air-Sea Interaction Regimes on July 1',30)
% Overlay the AS + BoB polygons used in panel (c)
hold on
roi_color = [0.30 0.30 0.30];
patch(AS_lon,  AS_lat,  'k','Edgecolor',roi_color, 'FaceColor','none', 'LineWidth',5,'linestyle','--');
patch(BoB_lon, BoB_lat, 'k','Edgecolor',roi_color, 'FaceColor','none', 'LineWidth',5,'linestyle','--');
hold off

ax4 = subplot(2,2,1);
pcolor(asi_regimes.lon,asi_regimes.lat,ws_clim_OAflux_short(:,:,dd)')
shading flat
makecb([0 15],cmocean('rain'),'Windspeed (m/s)',20,ax4);
hold on
geoshow(S, ...
    'Color', 'w', ...
    'LineWidth', 1)
hold off
makespruce(20)
ylabel('Latitude ($^\circ$N)', 'Interpreter','latex','fontsize',20)
xlabel('Longitude ($^\circ$E)', 'Interpreter','latex','fontsize',20)
maketitle('(a) 10m Windspeed Climatology on July 1',30)

% Overall title
sgtitle('', 'interpreter', 'latex', 'fontsize', 40);
annotation('textbox', [0.05 0.95 0.9 0.04], ...
    'String', 'Northern Indian Ocean ASI Regimes vs. Windspeed', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 40, 'interpreter', 'latex');

%% Spatial correlation between wind speed and ASI regimes

[C_ws, p_ws] = deal(nan(1,size(asi_regimes.param,3)));
for dd = 1:size(asi_regimes.param,3)
    ws_curr = ws_clim_OAflux_short(:,:,dd);
    reg_curr = asi_regimes.param(:,:,dd);
    goodind = ~isnan(ws_curr(:)) & ~isnan(reg_curr(:));
    [cmat,pmat] = corrcoef(ws_curr(goodind),reg_curr(goodind));
    C_ws(dd) = cmat(2,1);
    p_ws(dd) = pmat(2,1);
end

% Plot
subplot(2,2,4)
cplot = C_ws;
p_cutoff = BH(p_ws,0.1); % Benjamini-Hochberg procedure in time
cplot(p_ws>p_cutoff) = nan;
dateplot = datetime(2021,1,1):datetime(2021,12,31);  % non-leap year
plot(dateplot,cplot,'linewidth',2)
hold on
yline(0,'linewidth',2,'linestyle','--')
ylabel('Correlation','interpreter','latex')
makespruce(20)
maketitle('(d) Windspeed-ASI Regime Spatial Correlation',30)

% Adjust x-axis ticks and format
xticks(datetime(2021,1,1):calmonths(1):datetime(2021,12,31));
xtickformat('MMM');
ax = gca;
ax.XAxis.SecondaryLabel.Visible = 'off';  % the year is arbitrary for a climatology

fprintf('Maximum significant wind speed-ASI spatial correlation: R = %.2f on %s\n', ...
    max(cplot), datestr(dateplot(find(cplot==max(cplot),1)),'mmm dd'))

export_figure(fig1, P.fig_dir, 'asi_regimes_fig5_wind')
