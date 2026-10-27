% Nearshore vs. offshore ASI regimes in the Arabian Sea and Bay of Bengal
% (Section 3a: average P_o during March-May for grid points within two
% grid boxes of the coast vs. the interior)
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
% Outputs: printed MAM and MAMJJAS averages;
%          P.fig_dir/asi_regimes_nearshore_vs_offshore.png

clear; clc

P = asi_regimes_paths;

load(P.asi_OAFlux, 'asi_regimes')

%% Spatial mask: Arabian Sea + Bay of Bengal

% Build 2D lon/lat grid
[lon2d, lat2d] = meshgrid(asi_regimes.lon, asi_regimes.lat);
lon2d = lon2d';
lat2d = lat2d';

% Slanted AS boundary (5N,50E) to (26N,62E)
m = (62 - 50) / (26 - 5);
b = 50 - m * 5;

% Arabian Sea
AS_mask = ...
    lon2d >= 50 & lon2d <= 80 & ...
    lat2d >= 5  & lat2d <= 26 & ...
    lon2d >= (m * lat2d + b);

% Bay of Bengal
BoB_mask = ...
    lon2d >= 80 & lon2d <= 98 & ...
    lat2d >= 5  & lat2d <= 23;

% Final regional mask
region_mask = AS_mask | BoB_mask;

%% Distance to land

% Logical mask: 1 for land, 0 for ocean
land_logical = isnan(asi_regimes.param(:,:,1));

% Euclidean distance (in grid boxes) from each grid point to the nearest
% land point
[iland, jland] = find(land_logical);
[ii, jj] = ndgrid(1:size(land_logical,1), 1:size(land_logical,2));
dist_to_land = zeros(size(land_logical));
for k = 1:numel(ii)
    dist_to_land(k) = min(hypot(ii(k) - iland, jj(k) - jland));
end

% Define masks
nearshore_mask = dist_to_land <= 2;   % Within 2 grid boxes
offshore_mask = dist_to_land > 2;     % More than 2 grid boxes

%% Regional mean P_o by day of year

nearshore_meanparam = nan(1,365);
offshore_meanparam = nan(1,365);

for dd = 1:365
    asi_now = asi_regimes.param(:,:,dd);
    nearshore_vals = asi_now(region_mask&nearshore_mask);
    offshore_vals = asi_now(region_mask&offshore_mask);

    nearshore_meanparam(dd) = mean(nearshore_vals(:),'omitnan');
    offshore_meanparam(dd) = mean(offshore_vals(:),'omitnan');
end

%% Plot

% Day of year in the 365-day no-leap calendar; 2021 is a non-leap year
dates = datetime(2021,1,1:365);

fig1 = makefig(1800, 987);
hold on
plot(dates, nearshore_meanparam, 'b-', 'LineWidth', 3);
plot(dates, offshore_meanparam, 'r-', 'LineWidth', 3);
yline(0.5,'k--','LineWidth',2);
ylabel('ASI Regime Parameter','Interpreter','latex')
legend({'Nearshore','Offshore'},'Interpreter','latex','Location','best')
maketitle('Nearshore vs Offshore ASI Regimes',50)
makespruce(50)
hold off

export_figure(fig1, P.fig_dir, 'asi_regimes_nearshore_vs_offshore')

%% Seasonal averages

months = month(dates);

% Indices for March-May (MAM) and March-September (MAMJJAS)
mam_idx = months >= 3 & months <= 5;
mamjjas_idx = months >= 3 & months <= 9;

mam_nearshore_mean = mean(nearshore_meanparam(mam_idx),'omitnan');
mam_offshore_mean  = mean(offshore_meanparam(mam_idx),'omitnan');
mamjjas_nearshore_mean = mean(nearshore_meanparam(mamjjas_idx),'omitnan');
mamjjas_offshore_mean  = mean(offshore_meanparam(mamjjas_idx),'omitnan');

fprintf('MAM mean P_o:     nearshore %.2f, offshore %.2f\n', mam_nearshore_mean, mam_offshore_mean)
fprintf('MAMJJAS mean P_o: nearshore %.2f, offshore %.2f\n', mamjjas_nearshore_mean, mamjjas_offshore_mean)
