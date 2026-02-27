% Histograms of seasonal ASI regime parameter
% Alex Kinsella, February 2026

%% Load ASI regimes
load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_OAFlux_260210.mat')

param_all = asi_regimes.param;

% %% Spatial mask: Arabian Sea + Bay of Bengal
% 
% lonmatrix = repmat(asi_regimes.lon, [1 size(asi_regimes.param,[2 3])]);
% latmatrix = permute( ...
%     repmat(asi_regimes.lat, [1 size(asi_regimes.param,[1 3])]), ...
%     [2 1 3]);
% 
% % Northern Hemisphere mask
% NH_mask = latmatrix >= 0;
% 
% % Line separating Arabian Sea from land (5N,50E) to (26N,62E)
% m = (62 - 50) / (26 - 5);
% b = 50 - m * 5;
% 
% % Arabian Sea
% reg1_mask = ...
%     lonmatrix >= 50 & lonmatrix <= 80 & ...
%     latmatrix >= 0  & latmatrix <= 26 & ...
%     lonmatrix >= (m * latmatrix + b);
% 
% % Bay of Bengal
% reg2_mask = ...
%     lonmatrix >= 80 & lonmatrix <= 98 & ...
%     latmatrix >= 0  & latmatrix <= 23;
% 
% AS_BoB_mask = reg1_mask | reg2_mask;
% 
% final_mask = true(size(latmatrix));
% 
% if use_AS_BoB_mask
%     final_mask = final_mask & AS_BoB_mask;
% end
% 
% if use_NH_mask
%     final_mask = final_mask & NH_mask;
% end

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

%% Define bins
edges   = 0:0.05:1;
centers = edges(1:end-1) + diff(edges)/2;
bw      = diff(edges(1:2));

%% Season day-of-year indices
ind_FMAM = 32:151;       % Feb-Mar-Apr-May
ind_JJAS = 152:273;      % Jun-Jul-Aug-Sep
ind_ONDJ = [274:365 1:59];      % Oct-Nov-Dec-Jan

nt = size(param_all,3);

% Preallocate
param_masked = nan(size(param_all));

for dd = 1:nt
    slice = param_all(:,:,dd);
    slice(~region_mask) = NaN;
    param_masked(:,:,dd) = slice;
end

% Extract seasonal fields
param_FMAM = param_all(:,:,ind_FMAM);
param_JJAS = param_all(:,:,ind_JJAS);
param_ONDJ = param_all(:,:,ind_ONDJ);

param_FMAM = param_FMAM(~isnan(param_FMAM));
param_JJAS = param_JJAS(~isnan(param_JJAS));
param_ONDJ = param_ONDJ(~isnan(param_ONDJ));

%% Histogram counts 
c_FMAM = histcounts(param_FMAM(:), edges, 'Normalization','probability');
c_JJAS = histcounts(param_JJAS(:), edges, 'Normalization','probability');
c_ONDJ = histcounts(param_ONDJ(:), edges, 'Normalization','probability');

%% Figure
RGB = orderedcolors("gem");

% Side-by-side bar offsets within each bin
bar_width = 0.22;    % relative width (0–1)
offset    = bw/4;    % shift in x-units

fig = makefig;

hold on
b1 = bar(centers - offset, c_FMAM, bar_width, 'FaceColor', RGB(4,:), 'EdgeColor','none');
b2 = bar(centers, c_JJAS, bar_width, 'FaceColor', RGB(3,:), 'EdgeColor','none');
b3 = bar(centers + offset, c_ONDJ, bar_width, 'FaceColor', RGB(5,:), 'EdgeColor','none');
hold off

xlim([0 1])
ylim([0 0.25])
xlabel('ASI Parameter $P_o$','Interpreter','latex')
ylabel('Fractional Occurrence','Interpreter','latex')
legend([b1 b2 b3], {'Pre-Monsoon (FMAM)','Southwest Monsoon (JJAS)','Northeast Monsoon (ONDJ)'}, ...
    'Interpreter','latex','Location','northwest')

% Ocean / Atmosphere control annotations along Po (x) axis
text(0.98, -0.1, sprintf('%s\n%s','Ocean','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', ...
    'FontSize',25, ...
    'Interpreter','latex');

text(0.02, -0.1, sprintf('%s\n%s','Atm.','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',25, ...
    'Interpreter','latex');
makespruce(40);
maketitle('Seasonal Distributions of ASI Parameter $P_o$',50)
