% Figure 4: Histograms of seasonal ASI regime parameter
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% The histograms and percentages are computed over all grid points in the
% domain, as in the published Fig. 4 and Section 3a. (The published
% caption describes the Arabian Sea and Bay of Bengal regions only.)
%
% Inputs: P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
% Output: P.fig_dir/asi_regimes_fig4_histograms.png

clear; clc

P = asi_regimes_paths;

%% Load ASI regimes
load(P.asi_OAFlux, 'asi_regimes')

param_all = asi_regimes.param;

%% Define bins
edges   = 0:0.05:1;
centers = edges(1:end-1) + diff(edges)/2;
bw      = diff(edges(1:2));

%% Season day-of-year indices
ind_FMAM = 32:151;       % Feb-Mar-Apr-May
ind_JJAS = 152:273;      % Jun-Jul-Aug-Sep
ind_ONDJ = [274:365 1:31];      % Oct-Nov-Dec-Jan

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
bar_width = 0.22;    % relative width (0-1)
offset    = bw/4;    % shift in x-units

fig = makefig(1800, 987);  % window size of the published figure

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

%% Print results

fprintf(['FMAM: P_o<0.2 = %.1f%%, P_o>0.8 = %.1f%%\n' ...
'JJAS: P_o<0.2 = %.1f%%, P_o>0.8 = %.1f%%\n' ...
'ONDJ: P_o<0.2 = %.1f%%, P_o>0.8 = %.1f%%\n'], ...
100*mean(param_FMAM(:)<0.2), 100*mean(param_FMAM(:)>0.8), ...
100*mean(param_JJAS(:)<0.2), 100*mean(param_JJAS(:)>0.8), ...
100*mean(param_ONDJ(:)<0.2), 100*mean(param_ONDJ(:)>0.8));

export_figure(fig, P.fig_dir, 'asi_regimes_fig4_histograms')
