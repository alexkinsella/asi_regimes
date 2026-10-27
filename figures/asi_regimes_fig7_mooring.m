% Figure 7: Mooring heat budget analysis (Section 2c, Eqs. 9-14). Also
% produces supplemental Figure S1.
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  Arabian Sea Experiment mooring (P.AS_flux, P.AS_temp, P.AS_sal,
%          P.AS_curr) and ASIRI Bay of Bengal mooring (P.BoB_flux,
%          P.BoB_ts, P.BoB_curr); see README for sources
%          P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
%          P.oaflux_NIO (pipelines/asi_regimes_oaflux_processing.m), for
%          SST gradients
% Outputs: P.fig_dir/asi_regimes_fig7_mooring.png
%          P.fig_dir/asi_regimes_figS1_mooring_process.png
%
% Seawater density uses the GSW (TEOS-10) subset in util/gsw.

clear; clc

P = asi_regimes_paths;

%% ---------------- USER OPTIONS ----------------
f_asi = P.asi_OAFlux;

% Budget / ASI smoothing
smooth_days_budget = 7;
smooth_days_frac   = 30;
smooth_days_asi    = 30;

% Mechanism smoothing
smooth_days_horiz = 30;
smooth_days_vert  = 30;

% Unified MLD options
zref_mld      = 5;
drho_thresh   = 0.125;
min_mld       = 5;

% Minimum MLD deepening to count as finite entrainment event
min_deepen_m = 1;

% Climatological plot years
plot_year_common = 2015;

% BoB ASI comparison point
bob_alt_lon = 90;
bob_alt_lat = 15;

% ---------------- PLOT STYLE ----------------
colors.AS      = [0.000 0.447 0.741];  % blue
colors.BoB     = [0.850 0.325 0.098];  % orange
colors.ASI     = [0.10  0.10  0.10 ];  % near-black
colors.ASI_alt = [0.45  0.45  0.45 ];  % gray

lw_main = 2.3;
lw_aux  = 1.8;

% ---------------- FILE PATHS ----------------
% Arabian Sea
f_AS_flux = P.AS_flux;
f_AS_temp = P.AS_temp;
f_AS_sal  = P.AS_sal;
f_AS_curr = P.AS_curr;

% Bay of Bengal
f_BoB_flux = P.BoB_flux;
f_BoB_ts   = P.BoB_ts;
f_BoB_curr = P.BoB_curr;

% SST grid for gradients
f_sst_grid = P.oaflux_NIO;

% Velocity units
AS_vel_scale_to_mps  = 0.01;  % AS NetCDF says cm/s
BoB_vel_scale_to_mps = 1;     

% ========================================================================
% COMPUTE MOORING DIAGNOSTICS WITH MLD SENSITIVITY
% ========================================================================

% ASI box uncertainty
asi_box_half_width_deg = 1.5;  % 3 x 3 degree box

% ---------------- COMPUTE MOORING BUDGETS AND PROFILES ----------------
AS = calc_AS_all_v2(f_AS_flux, f_AS_temp, f_AS_sal, ...
                    zref_mld, drho_thresh, min_mld);

BoB = calc_BoB_all_v2(f_BoB_flux, f_BoB_ts, ...
                      zref_mld, drho_thresh, min_mld);

fprintf('\nUnified MLD criterion:\n')
fprintf('  Delta sigma = %.3f kg/m^3 relative to %.1f m\n', drho_thresh, zref_mld)
fprintf('  Minimum MLD = %.1f m\n', min_mld)
fprintf('  AS median MLD:  %.1f m\n', median(AS.MLD,'omitnan'))
fprintf('  BoB median MLD: %.1f m\n', median(BoB.MLD,'omitnan'))

% MLD lower/upper bounds from vertical density sampling uncertainty
AS.Hlo  = max(min_mld, AS.MLD - AS.MLD_unc);
AS.Hhi  = AS.MLD + AS.MLD_unc;
BoB.Hlo = max(min_mld, BoB.MLD - BoB.MLD_unc);
BoB.Hhi = BoB.MLD + BoB.MLD_unc;

% Preserve AS visual QC
AS.Hlo(AS.bad_mld_qc) = NaN;
AS.Hhi(AS.bad_mld_qc) = NaN;

% ---------------- MOORING OCEAN-GENERATED FRACTION ----------------
AS.frac      = calc_mooring_fraction_for_H(AS, AS.MLD, smooth_days_budget, smooth_days_frac);
AS.frac_Hlo  = calc_mooring_fraction_for_H(AS, AS.Hlo, smooth_days_budget, smooth_days_frac);
AS.frac_Hhi  = calc_mooring_fraction_for_H(AS, AS.Hhi, smooth_days_budget, smooth_days_frac);

BoB.frac     = calc_mooring_fraction_for_H(BoB, BoB.MLD, smooth_days_budget, smooth_days_frac);
BoB.frac_Hlo = calc_mooring_fraction_for_H(BoB, BoB.Hlo, smooth_days_budget, smooth_days_frac);
BoB.frac_Hhi = calc_mooring_fraction_for_H(BoB, BoB.Hhi, smooth_days_budget, smooth_days_frac);

AS.frac(AS.bad_mld_qc)     = NaN;
AS.frac_Hlo(AS.bad_mld_qc) = NaN;
AS.frac_Hhi(AS.bad_mld_qc) = NaN;

[AS.frac_lo, AS.frac_hi] = mld_envelope(AS.frac, AS.frac_Hlo, AS.frac_Hhi);

[BoB.frac_lo, BoB.frac_hi] = mld_envelope(BoB.frac, BoB.frac_Hlo, BoB.frac_Hhi);

% ---------------- LOAD / EXTRACT ASI WITH 3x3 DEGREE SPREAD ----------------
S = load(f_asi);
asi_regimes = S.asi_regimes;

% ---------------- HORIZONTAL MECHANISM: MIXED-LAYER MEAN SIGNED ADVECTION ----------------
AS.curr  = load_AS_current(f_AS_curr, AS_vel_scale_to_mps);
BoB.curr = load_BoB_current(f_BoB_curr, BoB_vel_scale_to_mps);

SST = load_flux_turb_NIO_sst(f_sst_grid);

AS.sstgrad  = extract_sst_gradient_at_point_flux_turb_NIO(SST, AS.lon, AS.lat, AS.date);
BoB.sstgrad = extract_sst_gradient_at_point_flux_turb_NIO(SST, BoB.lon, BoB.lat, BoB.date);

% Positive means cooling: u*dTdx + v*dTdy.
AS.horiz      = calc_mixed_layer_horizontal_cooling(AS.curr, AS, AS.MLD, smooth_days_horiz);
AS.horiz_Hlo  = calc_mixed_layer_horizontal_cooling(AS.curr, AS, AS.Hlo, smooth_days_horiz);
AS.horiz_Hhi  = calc_mixed_layer_horizontal_cooling(AS.curr, AS, AS.Hhi, smooth_days_horiz);

BoB.horiz     = calc_mixed_layer_horizontal_cooling(BoB.curr, BoB, BoB.MLD, smooth_days_horiz);
BoB.horiz_Hlo = calc_mixed_layer_horizontal_cooling(BoB.curr, BoB, BoB.Hlo, smooth_days_horiz);
BoB.horiz_Hhi = calc_mixed_layer_horizontal_cooling(BoB.curr, BoB, BoB.Hhi, smooth_days_horiz);

AS.horiz(AS.bad_mld_qc)     = NaN;
AS.horiz_Hlo(AS.bad_mld_qc) = NaN;
AS.horiz_Hhi(AS.bad_mld_qc) = NaN;

[AS.horiz_lo, AS.horiz_hi] = mld_envelope(AS.horiz, AS.horiz_Hlo, AS.horiz_Hhi);

[BoB.horiz_lo, BoB.horiz_hi] = mld_envelope(BoB.horiz, BoB.horiz_Hlo, BoB.horiz_Hhi);

fprintf('\nMedian mixed-layer horizontal advective cooling, positive = cooling:\n')
fprintf('  AS:  %.3g K/day\n', median(AS.horiz,'omitnan'))
fprintf('  BoB: %.3g K/day\n', median(BoB.horiz,'omitnan'))

% ---------------- VERTICAL ENTRAINMENT DIAGNOSTIC WITH MLD SENSITIVITY ----------------
AS.ent      = calc_finite_deepening_entrainment_for_H(AS, AS.MLD, min_deepen_m, smooth_days_vert);
AS.ent_Hlo  = calc_finite_deepening_entrainment_for_H(AS, AS.Hlo, min_deepen_m, smooth_days_vert);
AS.ent_Hhi  = calc_finite_deepening_entrainment_for_H(AS, AS.Hhi, min_deepen_m, smooth_days_vert);

BoB.ent     = calc_finite_deepening_entrainment_for_H(BoB, BoB.MLD, min_deepen_m, smooth_days_vert);
BoB.ent_Hlo = calc_finite_deepening_entrainment_for_H(BoB, BoB.Hlo, min_deepen_m, smooth_days_vert);
BoB.ent_Hhi = calc_finite_deepening_entrainment_for_H(BoB, BoB.Hhi, min_deepen_m, smooth_days_vert);

AS.ent(AS.bad_mld_qc)     = NaN;
AS.ent_Hlo(AS.bad_mld_qc) = NaN;
AS.ent_Hhi(AS.bad_mld_qc) = NaN;

[AS.ent_lo, AS.ent_hi] = mld_envelope(AS.ent, AS.ent_Hlo, AS.ent_Hhi);

[BoB.ent_lo, BoB.ent_hi] = mld_envelope(BoB.ent, BoB.ent_Hlo, BoB.ent_Hhi);

fprintf('\nMedian finite-deepening entrainment cooling, positive = cooling:\n')
fprintf('  AS:  %.3g K/day\n', median(AS.ent,'omitnan'))
fprintf('  BoB: %.3g K/day\n', median(BoB.ent,'omitnan'))

%% ========================================================================
%% PREPARE PLOTTING GRIDS
%% ========================================================================

% Common climatological display year
xgrid = (datetime(plot_year_common,1,1):datetime(plot_year_common,12,31))';
xnum  = datenum(xgrid);
nt = numel(xgrid);

month_ticks   = datetime(plot_year_common,1:12,1);
month_ticknum = datenum(month_ticks);
month_labs    = cellstr(datestr(month_ticks,'mmm'));

% Effective plotted record windows
AS_plot_keep  = AS.date  >= datetime(1994,10,3) & AS.date  < datetime(1995,9,15);
BoB_plot_keep = BoB.date >= datetime(2015,1,1)  & BoB.date < datetime(2016,1,1);

AS_plot_keep  = AS_plot_keep  & isfinite(AS.MLD);
BoB_plot_keep = BoB_plot_keep & isfinite(BoB.MLD);

doy_AS  = noleap_doy(AS.date);
doy_BoB = noleap_doy(BoB.date);

%% ---------------- Mooring ocean-fraction grids ----------------
AS_frac_grid    = grid_by_doy(AS.frac,     doy_AS,  AS_plot_keep, nt);
AS_frac_lo_grid = grid_by_doy(AS.frac_lo,  doy_AS,  AS_plot_keep, nt);
AS_frac_hi_grid = grid_by_doy(AS.frac_hi,  doy_AS,  AS_plot_keep, nt);

BoB_frac_grid    = grid_by_doy(BoB.frac,    doy_BoB, BoB_plot_keep, nt);
BoB_frac_lo_grid = grid_by_doy(BoB.frac_lo, doy_BoB, BoB_plot_keep, nt);
BoB_frac_hi_grid = grid_by_doy(BoB.frac_hi, doy_BoB, BoB_plot_keep, nt);

%% ---------------- Mechanism grids ----------------
AS_horiz_grid    = grid_by_doy(AS.horiz,    doy_AS,  AS_plot_keep, nt);
AS_horiz_lo_grid = grid_by_doy(AS.horiz_lo, doy_AS,  AS_plot_keep, nt);
AS_horiz_hi_grid = grid_by_doy(AS.horiz_hi, doy_AS,  AS_plot_keep, nt);

AS_ent_grid    = grid_by_doy(AS.ent,    doy_AS,  AS_plot_keep, nt);
AS_ent_lo_grid = grid_by_doy(AS.ent_lo, doy_AS,  AS_plot_keep, nt);
AS_ent_hi_grid = grid_by_doy(AS.ent_hi, doy_AS,  AS_plot_keep, nt);

BoB_horiz_grid    = grid_by_doy(BoB.horiz,    doy_BoB, BoB_plot_keep, nt);
BoB_horiz_lo_grid = grid_by_doy(BoB.horiz_lo, doy_BoB, BoB_plot_keep, nt);
BoB_horiz_hi_grid = grid_by_doy(BoB.horiz_hi, doy_BoB, BoB_plot_keep, nt);

BoB_ent_grid    = grid_by_doy(BoB.ent,    doy_BoB, BoB_plot_keep, nt);
BoB_ent_lo_grid = grid_by_doy(BoB.ent_lo, doy_BoB, BoB_plot_keep, nt);
BoB_ent_hi_grid = grid_by_doy(BoB.ent_hi, doy_BoB, BoB_plot_keep, nt);

%% ---------------- ASI annual-cycle grids ----------------
% Extract ASI directly on a canonical no-leap annual date vector to avoid
% duplicate DOY overwrites from mooring records spanning multiple years.
asi_plot_dates = datetime(plot_year_common,1,1) + days(0:nt-1);

asi_prctile_bounds = [10 90];

[AS_asi_grid, AS_asi_lo_grid, AS_asi_hi_grid] = extract_asi_box_prctile( ...
    asi_regimes, AS.lon, AS.lat, asi_plot_dates(:), ...
    asi_box_half_width_deg, asi_prctile_bounds);

[BoB_asi_grid, BoB_asi_lo_grid, BoB_asi_hi_grid] = extract_asi_box_prctile( ...
    asi_regimes, BoB.lon, BoB.lat, asi_plot_dates(:), ...
    asi_box_half_width_deg, asi_prctile_bounds);

[BoB_asi_alt_grid, BoB_asi_alt_lo_grid, BoB_asi_alt_hi_grid] = extract_asi_box_prctile( ...
    asi_regimes, bob_alt_lon, bob_alt_lat, asi_plot_dates(:), ...
    asi_box_half_width_deg, asi_prctile_bounds);

AS_asi_grid     = circular_movmean(AS_asi_grid,     smooth_days_asi);
AS_asi_lo_grid  = circular_movmean(AS_asi_lo_grid,  smooth_days_asi);
AS_asi_hi_grid  = circular_movmean(AS_asi_hi_grid,  smooth_days_asi);

BoB_asi_grid    = circular_movmean(BoB_asi_grid,    smooth_days_asi);
BoB_asi_lo_grid = circular_movmean(BoB_asi_lo_grid, smooth_days_asi);
BoB_asi_hi_grid = circular_movmean(BoB_asi_hi_grid, smooth_days_asi);

BoB_asi_alt_grid    = circular_movmean(BoB_asi_alt_grid,    smooth_days_asi);
BoB_asi_alt_lo_grid = circular_movmean(BoB_asi_alt_lo_grid, smooth_days_asi);
BoB_asi_alt_hi_grid = circular_movmean(BoB_asi_alt_hi_grid, smooth_days_asi);

AS_asi_lo_grid  = max(0, AS_asi_lo_grid);
AS_asi_hi_grid  = min(1, AS_asi_hi_grid);

BoB_asi_lo_grid = max(0, BoB_asi_lo_grid);
BoB_asi_hi_grid = min(1, BoB_asi_hi_grid);

BoB_asi_alt_lo_grid = max(0, BoB_asi_alt_lo_grid);
BoB_asi_alt_hi_grid = min(1, BoB_asi_alt_hi_grid);

%% ---------------- Figure style ----------------
cAS      = colors.AS;
cBoB     = colors.BoB;
cASI     = colors.ASI;
cASI_alt = colors.ASI_alt;

cHoriz   = [0.494 0.184 0.556];  % purple
cEnt     = [0.466 0.674 0.188];  % green

lw_moor = lw_main;
lw_asi  = lw_main;
lw_mech = lw_main;

fs_axis = 22;
fs_leg  = 25;

%% ========================================================================
%% FIGURE 7: MOORING HEAT BUDGET ANALYSIS
%% ========================================================================

fig7 = makefig(1800, 987);  % window size of the published figure

% ------------------------------------------------------------------------
% (a) Arabian Sea parameters
% ------------------------------------------------------------------------
ax1 = subplot(1,2,1);

shade_between(xnum, AS_frac_lo_grid, AS_frac_hi_grid, cAS, 0.16); 
hold on
shade_between(xnum, AS_asi_lo_grid, AS_asi_hi_grid, cASI, 0.12);

hAS = plot(xnum, AS_frac_grid, ...
           'Color', cAS, ...
           'LineWidth', lw_moor);

hPo = plot(xnum, AS_asi_grid, ...
           'Color', cASI, ...
           'LineWidth', lw_asi);

ylim([0 1])
xlim([xnum(1) xnum(end)])
xticks(month_ticknum)
xticklabels(month_labs)

makeylabel('Parameter Value ($P_o$ or $P_\mathrm{moor}$)')
grid on
box on
set(ax1, 'FontSize', fs_axis, 'TickDir','out', 'LineWidth',1.2)

makespruce(25)
maketitle('(a) Arabian Sea', 30)

% ------------------------------------------------------------------------
% (b) Bay of Bengal parameters
% ------------------------------------------------------------------------
ax2 = subplot(1,2,2);

shade_between(xnum, BoB_frac_lo_grid, BoB_frac_hi_grid, cBoB, 0.16); 
hold on
shade_between(xnum, BoB_asi_lo_grid, BoB_asi_hi_grid, cASI, 0.12);
%shade_between(xnum, BoB_asi_alt_lo_grid, BoB_asi_alt_hi_grid, cASI_alt, 0.08);

hBoB = plot(xnum, BoB_frac_grid, ...
            'Color', cBoB, ...
            'LineWidth', lw_moor);

plot(xnum, BoB_asi_grid, ...
     'Color', cASI, ...
     'LineWidth', lw_asi);

hPoAlt = plot(xnum, BoB_asi_alt_grid, '--', ...
              'Color', cASI_alt, ...
              'LineWidth', lw_aux);

ylim([0 1])
xlim([xnum(1) xnum(end)])
xticks(month_ticknum)
xticklabels(month_labs)

%makeylabel('$P_o$ or $P_\mathrm{moor}$')
grid on
box on
set(ax2, 'FontSize', fs_axis, 'TickDir','out', 'LineWidth',1.2)

makespruce(25)
maketitle('(b) Bay of Bengal', 30)

% ---- make room for bottom shared legend ----
% Raise the bottom of each axis, but keep the top nearly fixed.
raise_bottom = 0.075;   % creates room below
lower_top    = 0.075;   % creates a little extra room below sgtitle

p1 = ax1.Position;
p2 = ax2.Position;

new_y1 = p1(2) + raise_bottom;
new_y2 = p2(2) + raise_bottom;

new_h1 = p1(4) - raise_bottom - lower_top;
new_h2 = p2(4) - raise_bottom - lower_top;

ax1.Position = [p1(1), new_y1, p1(3), new_h1];
ax2.Position = [p2(1), new_y2, p2(3), new_h2];

% ---- Figure 7 legend ----
legax7 = axes(fig7, ...
              'Position', [0.08 0.005 0.84 0.045], ...
              'Visible', 'off');
hold(legax7, 'on')

l1 = plot(legax7, nan, nan, '-',  'Color', cAS,      'LineWidth', lw_moor);
l2 = plot(legax7, nan, nan, '-',  'Color', cBoB,     'LineWidth', lw_moor);
l3 = plot(legax7, nan, nan, '-',  'Color', cASI,     'LineWidth', lw_asi);
l4 = plot(legax7, nan, nan, '--', 'Color', cASI_alt, 'LineWidth', lw_aux);

legend(legax7, [l1 l2 l3 l4], ...
       {'AS mooring $P_\mathrm{moor}$', ...
        'BoB mooring $P_\mathrm{moor}$', ...
        'ASI parameter $P_o$ at mooring', ...
        '$P_o$ at 15$^\circ$N, 90$^\circ$E'}, ...
       'Orientation','horizontal', ...
       'Location','south', ...
       'Interpreter','latex', ...
       'FontSize',fs_leg, ...
       'Box','on')

sgtitle('Mooring SST Tendency and ASI Regime Comparison', ...
        'Interpreter','latex', ...
        'FontSize',35)

export_figure(fig7, P.fig_dir, 'asi_regimes_fig7_mooring')

%% ========================================================================
%% FIGURE S1: MOORING OCEAN-PROCESS MECHANISMS
%% ========================================================================

% ---- convert mechanism diagnostics from cooling-positive to warming-positive ----
% Current convention:
%   AS_horiz_grid, AS_ent_grid, etc. are positive for cooling.
%
% Plotting convention for Fig. S1:
%   positive = warming, negative = cooling.

AS_horiz_warm_grid = -AS_horiz_grid;
BoB_horiz_warm_grid = -BoB_horiz_grid;

AS_ent_warm_grid = -AS_ent_grid;
BoB_ent_warm_grid = -BoB_ent_grid;

% Envelopes: sign flip swaps lower/upper bounds.
AS_horiz_warm_lo_grid = -AS_horiz_hi_grid;
AS_horiz_warm_hi_grid = -AS_horiz_lo_grid;

BoB_horiz_warm_lo_grid = -BoB_horiz_hi_grid;
BoB_horiz_warm_hi_grid = -BoB_horiz_lo_grid;

AS_ent_warm_lo_grid = -AS_ent_hi_grid;
AS_ent_warm_hi_grid = -AS_ent_lo_grid;

BoB_ent_warm_lo_grid = -BoB_ent_hi_grid;
BoB_ent_warm_hi_grid = -BoB_ent_lo_grid;

figS1 = makefig(1800, 987);  % window size of the published figure

% ------------------------------------------------------------------------
% (a) Arabian Sea mechanisms
% ------------------------------------------------------------------------
axS1 = subplot(1,2,1);

hold on

% Horizontal-advection shading intentionally omitted because its dominant
% uncertainty is unresolved SST-gradient structure, not MLD uncertainty.
% shade_between(xnum, AS_horiz_warm_lo_grid, AS_horiz_warm_hi_grid, cHoriz, 0.16);

shade_between(xnum, AS_ent_warm_lo_grid, AS_ent_warm_hi_grid, cEnt, 0.16);

hH = plot(xnum, AS_horiz_warm_grid, '-', ...
          'Color', cHoriz, ...
          'LineWidth', lw_mech);

hE = plot(xnum, AS_ent_warm_grid, '-', ...
          'Color', cEnt, ...
          'LineWidth', lw_mech);

yline(0, 'k-', 'HandleVisibility','off')

xlim([xnum(1) xnum(end)])
xticks(month_ticknum)
xticklabels(month_labs)

legend([hH hE], ...
       {'Horizontal advection', 'Vertical entrainment'}, ...
       'Interpreter','latex', ...
       'Location','best')

makeylabel('Mixed layer temperature tendency (K day$^{-1}$)')
grid on
box on
set(axS1, 'FontSize', fs_axis, 'TickDir','out', 'LineWidth',1.2)

makespruce(20)
maketitle('(a) Arabian Sea', 30)

% ------------------------------------------------------------------------
% (b) Bay of Bengal mechanisms
% ------------------------------------------------------------------------
axS2 = subplot(1,2,2);

hold on

% Horizontal-advection shading intentionally omitted because its dominant
% uncertainty is unresolved SST-gradient structure, not MLD uncertainty.
% shade_between(xnum, BoB_horiz_warm_lo_grid, BoB_horiz_warm_hi_grid, cHoriz, 0.16);

shade_between(xnum, BoB_ent_warm_lo_grid, BoB_ent_warm_hi_grid, cEnt, 0.16);

hH2 = plot(xnum, BoB_horiz_warm_grid, '-', ...
           'Color', cHoriz, ...
           'LineWidth', lw_mech);

hE2 = plot(xnum, BoB_ent_warm_grid, '-', ...
           'Color', cEnt, ...
           'LineWidth', lw_mech);

yline(0, 'k-', 'HandleVisibility','off')

xlim([xnum(1) xnum(end)])
xticks(month_ticknum)
xticklabels(month_labs)

makeylabel('Mixed layer temperature tendency (K day$^{-1}$)')
grid on
box on
set(axS2, 'FontSize', fs_axis, 'TickDir','out', 'LineWidth',1.2)

makespruce(20)
maketitle('(b) Bay of Bengal', 30)

sgtitle('Mooring-Based Ocean Mixed Layer Process Diagnostics', ...
        'Interpreter','latex', ...
        'FontSize',40)

export_figure(figS1, P.fig_dir, 'asi_regimes_figS1_mooring_process')

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function C = load_AS_current(fname, scale_to_mps)

    C.name = 'Arabian Sea currents';
    C.t = epic_time(fname);
    C.z = double(ncread(fname,'depth'));

    C.u = double(squeeze(ncread(fname,'east')));
    C.v = double(squeeze(ncread(fname,'north')));

    C.t = C.t(:);
    C.z = C.z(:);

    if size(C.u,1) ~= numel(C.t); C.u = C.u'; end
    if size(C.v,1) ~= numel(C.t); C.v = C.v'; end

    C.u(abs(C.u) > 1e20) = NaN;
    C.v(abs(C.v) > 1e20) = NaN;

    C.u = C.u * scale_to_mps;
    C.v = C.v * scale_to_mps;
end

function C = load_BoB_current(fname, scale_to_mps)

    S = load(fname);

    C.name = 'Bay of Bengal currents';
    C.t = get_bob_time(S.mday, S.yday, S.year);
    C.z = double(S.depth(:));

    C.u = double(S.east);
    C.v = double(S.north);

    if size(C.u,1) ~= numel(C.t); C.u = C.u'; end
    if size(C.v,1) ~= numel(C.t); C.v = C.v'; end

    C.u(abs(C.u) > 1e20) = NaN;
    C.v(abs(C.v) > 1e20) = NaN;

    C.u = C.u * scale_to_mps;
    C.v = C.v * scale_to_mps;
end

function SST = load_flux_turb_NIO_sst(fname)

    S = load(fname, 'flux_turb_NIO');
    F = S.flux_turb_NIO;

    SST.lon = double(F.lon(:));
    SST.lat = double(F.lat(:));
    SST.t   = double(F.time(:));
    SST.sst = double(F.sst);

    SST.sst(abs(SST.sst) > 1e20) = NaN;

end

function G = extract_sst_gradient_at_point_flux_turb_NIO(SST, lon0, lat0, dates)

    lon = SST.lon(:);
    lat = SST.lat(:);
    t_sst = SST.t(:);
    sst = SST.sst;

    lon0_use = lon0;
    if max(lon) > 180 && lon0 < 0
        lon0_use = mod(lon0, 360);
    elseif max(lon) <= 180 && lon0 > 180
        lon0_use = lon0 - 360;
    end

    [~, ilon] = min(abs(lon - lon0_use));
    [~, ilat] = min(abs(lat - lat0));

    t_target = datenum(dates(:));

    Tw_ts = squeeze(sst(ilon-1, ilat,   :));
    Te_ts = squeeze(sst(ilon+1, ilat,   :));
    Ts_ts = squeeze(sst(ilon,   ilat-1, :));
    Tn_ts = squeeze(sst(ilon,   ilat+1, :));

    Re = 6371000;
    dlon_rad = deg2rad(lon(ilon+1) - lon(ilon-1));
    dlat_rad = deg2rad(lat(ilat+1) - lat(ilat-1));

    dx = Re * cosd(lat(ilat)) * dlon_rad;
    dy = Re * dlat_rad;

    dTdx_ts = (Te_ts - Tw_ts) ./ dx;
    dTdy_ts = (Tn_ts - Ts_ts) ./ dy;

    G.dTdx = interp1(t_sst, dTdx_ts, t_target, 'linear', NaN);
    G.dTdy = interp1(t_sst, dTdy_ts, t_target, 'linear', NaN);
    G.gradmag = sqrt(G.dTdx.^2 + G.dTdy.^2);

    G.ilon = ilon;
    G.ilat = ilat;
    G.lon_grid = lon(ilon);
    G.lat_grid = lat(ilat);
end

function Tbar = layer_mean_between_depths(Tprof, z, z1, z2)

    Tbar = NaN;

    if ~isfinite(z1) || ~isfinite(z2) || z2 <= z1
        return
    end

    z = z(:);
    Tprof = Tprof(:);

    good = isfinite(z) & isfinite(Tprof);
    if sum(good) < 2
        return
    end

    zz = z(good);
    TT = Tprof(good);

    [zz,I] = sort(zz);
    TT = TT(I);

    if z1 < min(zz) || z2 > max(zz)
        return
    end

    T1 = interp1(zz, TT, z1, 'linear');
    T2 = interp1(zz, TT, z2, 'linear');

    use = zz > z1 & zz < z2;

    z_int = [z1; zz(use); z2];
    T_int = [T1; TT(use); T2];

    Tbar = trapz(z_int, T_int) / (z2 - z1);
end

function T_on_newz = interp_profiles_to_depths(T, z_old, z_new, max_extrap)

    z_old = z_old(:);
    z_new = z_new(:);

    T_on_newz = nan(size(T,1), numel(z_new));

    for ii = 1:size(T,1)
        Ti = T(ii,:).';
        good = isfinite(Ti) & isfinite(z_old);

        if sum(good) < 2
            continue
        end

        zg = z_old(good);
        Tg = Ti(good);

        [zg,I] = sort(zg);
        Tg = Tg(I);

        Ttmp = interp1(zg, Tg, z_new, 'linear', NaN);

        shallow = isnan(Ttmp) & z_new < min(zg) & (min(zg)-z_new) <= max_extrap;
        deep    = isnan(Ttmp) & z_new > max(zg) & (z_new-max(zg)) <= max_extrap;

        Ttmp(shallow) = Tg(1);
        Ttmp(deep)    = Tg(end);

        T_on_newz(ii,:) = Ttmp;
    end
end

function Tml = calc_Tml_from_profile(T, z, MLD, Tsurf, zsurf)

    z = z(:);

    if nargin < 4 || isempty(Tsurf)
        Tsurf = nan(size(MLD));
    end

    if nargin < 5 || isempty(zsurf)
        zsurf = 0;
    end

    Tml = nan(numel(MLD),1);

    for ii = 1:numel(MLD)
        h = MLD(ii);

        if ~isfinite(h) || h <= 0
            continue
        end

        Ti = T(ii,:).';
        good = isfinite(Ti) & isfinite(z);

        if sum(good) < 2
            continue
        end

        zg = z(good);
        Tg = Ti(good);

        [zg,I] = sort(zg);
        Tg = Tg(I);

        z_prof = zg(:);
        T_prof = Tg(:);

        if ~isnan(Tsurf(ii))
            z_prof = [zsurf; z_prof];
            T_prof = [Tsurf(ii); T_prof];
        end

        [z_prof,I] = sort(z_prof);
        T_prof = T_prof(I);

        [z_prof,ia] = unique(z_prof);
        T_prof = T_prof(ia);

        if min(z_prof) > 0
            z_prof = [0; z_prof];
            T_prof = [T_prof(1); T_prof];
        end

        if h <= max(z_prof)
            T_h = interp1(z_prof, T_prof, h, 'linear');
        else
            T_h = T_prof(end);
        end

        z_int = [z_prof; h];
        T_int = [T_prof; T_h];

        use = z_int >= 0 & z_int <= h & isfinite(T_int);

        if sum(use) < 2
            continue
        end

        [zsort,I] = sort(z_int(use));
        Tsort = T_int(use);
        Tsort = Tsort(I);

        Tml(ii) = trapz(zsort, Tsort) / h;
    end
end

function doy = noleap_doy(dates)

    doy = day(dates,'dayofyear');
    yr = year(dates);

    is_leap = eomday(yr,2) == 29;
    after_feb28 = dates > datetime(yr,2,28);

    is_feb29 = month(dates) == 2 & day(dates) == 29;
    doy(is_feb29) = NaN;

    shift = is_leap & after_feb28 & ~is_feb29;
    doy(shift) = doy(shift) - 1;
end

function ysm = circular_movmean(y, win)

    y = y(:);
    n = numel(y);

    if n == 0
        ysm = y;
        return
    end

    ypad = [y; y; y];
    ypad_sm = movmean(ypad, win, 'omitnan');
    ysm = ypad_sm(n+1:2*n);
end

function t = epic_time(fname)

    time  = double(ncread(fname,'time'));
    time2 = double(ncread(fname,'time2'));
    t = time - 1721059 + time2/(1000*86400);
    t = t(:);
end

function t = get_bob_time(mday, yday, year0)

    mday = mday(:);
    yday = yday(:);

    if median(mday,'omitnan') > 1e5
        t = double(mday);
    else
        t = datenum(year0,1,0) + double(yday);
    end
end

function Xd = daily_bin(t, X, td)

    t = t(:);

    if isvector(X)
        X = X(:);
    end

    Xd = nan(numel(td), size(X,2));

    for ii = 1:numel(td)
        ind = t >= td(ii)-0.5 & t < td(ii)+0.5;
        if any(ind)
            Xd(ii,:) = mean(X(ind,:), 1, 'omitnan');
        end
    end
end

function out = calc_AS_all_v2(f_flux, f_temp, f_sal, zref, drho_thresh, min_mld)

    tT = epic_time(f_temp);
    tS = epic_time(f_sal);
    tQ = epic_time(f_flux);

    zT = double(ncread(f_temp,'depth'));
    zS = double(ncread(f_sal,'depth'));

    lon = double(ncread(f_temp,'lon'));
    lat = double(ncread(f_temp,'lat'));

    T = double(squeeze(ncread(f_temp,'temp')));
    S = double(squeeze(ncread(f_sal,'sal')));

    Qnet = double(squeeze(ncread(f_flux,'QN')));
    QH   = double(squeeze(ncread(f_flux,'QH')));
    QB   = double(squeeze(ncread(f_flux,'QB')));
    Qs   = double(squeeze(ncread(f_flux,'Qs')));

    if size(T,1) ~= numel(tT); T = T'; end
    if size(S,1) ~= numel(tS); S = S'; end

    T(abs(T) > 1e20) = NaN;
    S(abs(S) > 1e20) = NaN;
    Qnet(abs(Qnet) > 1e20) = NaN;
    QH(abs(QH) > 1e20) = NaN;
    QB(abs(QB) > 1e20) = NaN;
    Qs(abs(Qs) > 1e20) = NaN;

    td = (ceil(max([min(tT), min(tS), min(tQ)])):floor(min([max(tT), max(tS), max(tQ)])))';
    date_d = datetime(td,'ConvertFrom','datenum');

    T_d = daily_bin(tT, T, td);
    S_d = daily_bin(tS, S, td);

    Qnet_d = daily_bin(tQ, Qnet, td);
    QH_d   = daily_bin(tQ, QH,   td);
    QB_d   = daily_bin(tQ, QB,   td);
    Qs_d   = daily_bin(tQ, Qs,   td);

    T_on_S = interp_profiles_to_depths(T_d, zT, zS, 10);

    % TEOS-10 potential density (reference pressure 0); depth in m is
    % used as pressure in dbar
    SA = gsw_SA_from_SP(S_d, zS(:)', lon, lat);
    CT = gsw_CT_from_t(SA, T_on_S, zS(:)');
    rho = gsw_rho(SA, CT, 0);

    [MLD, MLD_unc] = calc_mld_density_ref_unc(rho, zS, zref, drho_thresh);
    MLD = max(MLD, min_mld);

    bad_mld_qc = (td >= 728766 & td <= 728770) | (td >= 728912);
    MLD(bad_mld_qc) = NaN;
    MLD_unc(bad_mld_qc) = NaN;

    out.name = 'Arabian Sea';
    out.lon = lon;
    out.lat = lat;
    out.td = td;
    out.date = date_d;
    out.zT = zT(:);
    out.zS = zS(:);
    out.T_d = T_d;
    out.S_d = S_d;
    out.rho_d = rho;
    out.MLD = MLD;
    out.MLD_unc = MLD_unc;
    out.Qnet_d = Qnet_d;
    out.Qturb_d = QH_d + QB_d;
    out.Qs_d = Qs_d;
    out.Tsurf_d = nan(size(MLD));
    out.zsurf = [];
    out.bad_mld_qc = bad_mld_qc;

    out.sstgrad = [];
end

function out = calc_BoB_all_v2(f_flux, f_ts, zref, drho_thresh, min_mld)

    flux = load(f_flux);
    ts   = load(f_ts);

    t_flux = get_bob_time(flux.mday, flux.yday, flux.year);
    t_ts   = get_bob_time(ts.mday,   ts.yday,   ts.year);

    Qnet = flux.QN(:);
    QH   = flux.QH(:);
    QB   = flux.QB(:);
    Qs   = flux.Qs(:);
    tskin = flux.tskin(:);

    z = double(ts.depth(:));
    T = double(ts.temp);
    sig = double(ts.sigma);

    if size(T,1) ~= numel(t_ts); T = T'; end
    if size(sig,1) ~= numel(t_ts); sig = sig'; end

    Qnet(abs(Qnet) > 1e20) = NaN;
    QH(abs(QH) > 1e20) = NaN;
    QB(abs(QB) > 1e20) = NaN;
    Qs(abs(Qs) > 1e20) = NaN;
    tskin(abs(tskin) > 1e20) = NaN;
    T(abs(T) > 1e20) = NaN;
    sig(abs(sig) > 1e20) = NaN;

    t = t_ts;

    Qnet_i  = interp1(t_flux, Qnet,  t, 'linear', NaN);
    QH_i    = interp1(t_flux, QH,    t, 'linear', NaN);
    QB_i    = interp1(t_flux, QB,    t, 'linear', NaN);
    Qs_i    = interp1(t_flux, Qs,    t, 'linear', NaN);
    tskin_i = interp1(t_flux, tskin, t, 'linear', NaN);

    td = (ceil(min(t)):floor(max(t)))';
    date_d = datetime(td,'ConvertFrom','datenum');

    T_d   = daily_bin(t, T, td);
    sig_d = daily_bin(t, sig, td);

    Qnet_d  = daily_bin(t, Qnet_i,  td);
    QH_d    = daily_bin(t, QH_i,    td);
    QB_d    = daily_bin(t, QB_i,    td);
    Qs_d    = daily_bin(t, Qs_i,    td);
    tskin_d = daily_bin(t, tskin_i, td);

    [MLD, MLD_unc] = calc_mld_density_ref_unc(sig_d, z, zref, drho_thresh);
    MLD = max(MLD, min_mld);

    if isfield(flux,'zsst_in')
        zsst = flux.zsst_in;
    else
        zsst = 0;
    end

    out.name = 'Bay of Bengal';
    out.lon = ts.lon;
    out.lat = ts.lat;
    out.td = td;
    out.date = date_d;
    out.zT = z;
    out.zS = z;
    out.T_d = T_d;
    out.rho_d = sig_d + 1000;
    out.MLD = MLD;
    out.MLD_unc = MLD_unc;
    out.Qnet_d = Qnet_d;
    out.Qturb_d = QH_d + QB_d;
    out.Qs_d = Qs_d;
    out.Tsurf_d = tskin_d;
    out.zsurf = zsst;
    out.bad_mld_qc = false(size(MLD));

    out.sstgrad = [];
end

function frac = calc_mooring_fraction_for_H(M, H, smooth_days_budget, smooth_days_frac)

    rho0 = 1025;
    cp   = 3990;

    sw_R     = 0.58;
    sw_zeta1 = 0.35;
    sw_zeta2 = 23.0;

    H = H(:);

    Tml = calc_Tml_from_profile(M.T_d, M.zT, H, M.Tsurf_d, M.zsurf);

    dTdt_obs = nan(size(Tml));
    td = M.td(:);

    dTdt_obs(2:end-1) = (Tml(3:end) - Tml(1:end-2)) ./ ...
                        ((td(3:end) - td(1:end-2))*86400);

    frac_pen = sw_R .* exp(-H ./ sw_zeta1) + ...
               (1 - sw_R) .* exp(-H ./ sw_zeta2);

    Qpen_d = M.Qs_d .* frac_pen;
    Qnet_abs_d = M.Qnet_d - Qpen_d;

    dTdt_turb  = M.Qturb_d    ./ (rho0 * cp .* H);
    dTdt_surf  = Qnet_abs_d   ./ (rho0 * cp .* H);
    dTdt_ocean = dTdt_obs - dTdt_surf;

    turb_Kday  = dTdt_turb  * 86400;
    ocean_Kday = dTdt_ocean * 86400;

    turb_sm  = movmean(turb_Kday,  smooth_days_budget, 'omitnan');
    ocean_sm = movmean(ocean_Kday, smooth_days_budget, 'omitnan');

    frac = abs(ocean_sm) ./ sqrt(turb_sm.^2 + ocean_sm.^2);
    frac = movmean(frac, smooth_days_frac, 'omitnan');

    frac(~isfinite(H)) = NaN;
    frac(M.bad_mld_qc) = NaN;
end

function horiz_sm = calc_mixed_layer_horizontal_cooling(C, M, H, smooth_days)

    u_d = daily_bin(C.t, C.u, M.td);
    v_d = daily_bin(C.t, C.v, M.td);

    H = H(:);
    n = numel(H);

    u_ml = nan(n,1);
    v_ml = nan(n,1);

    for ii = 1:n
        if ~isfinite(H(ii)) || H(ii) <= 0
            continue
        end

        u_ml(ii) = layer_mean_profile_extrap(C.z, u_d(ii,:), H(ii));
        v_ml(ii) = layer_mean_profile_extrap(C.z, v_d(ii,:), H(ii));
    end

    % Positive = cooling because dT/dt_adv = -(u dTdx + v dTdy).
    horiz_cooling_Kday = ...
        (u_ml .* M.sstgrad.dTdx(:) + v_ml .* M.sstgrad.dTdy(:)) * 86400;

    horiz_sm = movmean(horiz_cooling_Kday, smooth_days, 'omitnan');
    horiz_sm(~isfinite(H)) = NaN;
    horiz_sm(M.bad_mld_qc) = NaN;
end

function ent_sm = calc_finite_deepening_entrainment_for_H(M, H, min_deepen_m, smooth_days)

    H = H(:);
    Tml = calc_Tml_from_profile(M.T_d, M.zT, H, M.Tsurf_d, M.zsurf);
    t = datenum(M.date(:));

    n = numel(H);
    cooling_Kday = nan(n,1);

    for ii = 2:n

        if any(~isfinite([H(ii-1), H(ii), Tml(ii-1), t(ii-1), t(ii)]))
            continue
        end

        dt_day = t(ii) - t(ii-1);
        if dt_day <= 0
            continue
        end

        dH = H(ii) - H(ii-1);

        if dH <= min_deepen_m
            cooling_Kday(ii) = 0;
            continue
        end

        H0 = H(ii-1);
        H1 = H(ii);

        T_layer = layer_mean_between_depths(M.T_d(ii-1,:), M.zT, H0, H1);

        if ~isfinite(T_layer)
            continue
        end

        T_new = (H0*Tml(ii-1) + (H1-H0)*T_layer) / H1;

        cooling_Kday(ii) = (Tml(ii-1) - T_new) / dt_day;
    end

    ent_sm = movmean(cooling_Kday, smooth_days, 'omitnan');
    ent_sm(~isfinite(H)) = NaN;
    ent_sm(M.bad_mld_qc) = NaN;
end

function [MLD, MLD_unc] = calc_mld_density_ref_unc(rho_or_sigma, z, zref, drho_thresh)

    z = z(:);
    MLD = nan(size(rho_or_sigma,1),1);
    MLD_unc = nan(size(rho_or_sigma,1),1);

    max_ref_extrap = 3;
    min_points = 3;

    for ii = 1:size(rho_or_sigma,1)

        r = rho_or_sigma(ii,:).';
        good = isfinite(r) & isfinite(z);

        if sum(good) < min_points
            continue
        end

        zz = z(good);
        rr = r(good);

        [zz,I] = sort(zz);
        rr = rr(I);

        if zref < min(zz)
            if min(zz)-zref <= max_ref_extrap
                rref = rr(1);
            else
                continue
            end
        elseif zref > max(zz)
            if zref-max(zz) <= max_ref_extrap
                rref = rr(end);
            else
                continue
            end
        else
            rref = interp1(zz, rr, zref, 'linear');
        end

        below = zz >= zref;
        zzb = zz(below);
        rrb = rr(below);

        if numel(zzb) < 2
            continue
        end

        drho = rrb - rref;
        kk = find(drho >= drho_thresh, 1, 'first');

        if isempty(kk)
            MLD(ii) = max(zzb);
            MLD_unc(ii) = 0.5 * abs(zzb(end) - zzb(end-1));

        elseif kk == 1
            MLD(ii) = zzb(1);
            MLD_unc(ii) = 0.5 * abs(zzb(2) - zzb(1));

        else
            z1 = zzb(kk-1); z2 = zzb(kk);
            d1 = drho(kk-1); d2 = drho(kk);

            if d2 ~= d1
                MLD(ii) = z1 + (drho_thresh-d1)*(z2-z1)/(d2-d1);
            else
                MLD(ii) = z2;
            end

            MLD_unc(ii) = 0.5 * abs(z2 - z1);
        end
    end
end

function ubar = layer_mean_profile_extrap(z, u, H)

    z = z(:);
    u = u(:);

    good = isfinite(z) & isfinite(u);
    if sum(good) < 1 || ~isfinite(H) || H <= 0
        ubar = NaN;
        return
    end

    zz = z(good);
    uu = u(good);

    [zz,I] = sort(zz);
    uu = uu(I);

    % Need a profile from 0 to H.
    z_inside = zz(zz > 0 & zz < H);
    z_int = unique([0; z_inside; H]);

    u_int = nan(size(z_int));

    for jj = 1:numel(z_int)
        zq = z_int(jj);

        if zq <= min(zz)
            % assume shallowest measured current extends to surface
            u_int(jj) = uu(1);

        elseif zq >= max(zz)
            % extrapolate below deepest available velocity using shear at base
            if numel(zz) >= 2
                shear = (uu(end) - uu(end-1)) / (zz(end) - zz(end-1));
                u_int(jj) = uu(end) + shear * (zq - zz(end));
            else
                u_int(jj) = uu(end);
            end

        else
            u_int(jj) = interp1(zz, uu, zq, 'linear');
        end
    end

    ubar = trapz(z_int, u_int) / H;
end

function [asi_center, asi_lo, asi_hi] = extract_asi_box_prctile(asi_regimes, lon0, lat0, dates, half_width_deg, prctile_bounds)

    if nargin < 6 || isempty(prctile_bounds)
        prctile_bounds = [10 90];
    end

    param = asi_regimes.param;

    lon = [];
    lat = [];

    for name = {'lon','longitude','lons'}
        if isfield(asi_regimes, name{1})
            lon = asi_regimes.(name{1});
            break
        end
    end

    for name = {'lat','latitude','lats'}
        if isfield(asi_regimes, name{1})
            lat = asi_regimes.(name{1});
            break
        end
    end

    lon = lon(:);
    lat = lat(:);

    lon0_use = lon0;
    if max(lon) > 180 && lon0 < 0
        lon0_use = mod(lon0, 360);
    elseif max(lon) <= 180 && lon0 > 180
        lon0_use = lon0 - 360;
    end

    P = orient_asi_to_lon_lat_time(param, lon, lat);

    [~, ilon0] = min(abs(lon - lon0_use));
    [~, ilat0] = min(abs(lat - lat0));

    ilon_box = find(abs(lon - lon0_use) <= half_width_deg);
    ilat_box = find(abs(lat - lat0)     <= half_width_deg);

    doy = noleap_doy(dates);

    asi_center = nan(size(doy));
    asi_lo     = nan(size(doy));
    asi_hi     = nan(size(doy));

    for ii = 1:numel(doy)
        if isnan(doy(ii)) || doy(ii) < 1 || doy(ii) > size(P,3)
            continue
        end

        slice = P(:,:,doy(ii));

        asi_center(ii) = slice(ilon0, ilat0);

        boxvals = slice(ilon_box, ilat_box);
        boxvals = boxvals(:);
        boxvals = boxvals(isfinite(boxvals));

        if ~isempty(boxvals)
            asi_lo(ii) = prctile(boxvals, prctile_bounds(1));
            asi_hi(ii) = prctile(boxvals, prctile_bounds(2));
        end
    end
end

% function [asi_center, asi_std] = extract_asi_box_stats(asi_regimes, lon0, lat0, dates, half_width_deg)
% 
%     param = asi_regimes.param;
% 
%     lon = [];
%     lat = [];
% 
%     for name = {'lon','longitude','lons'}
%         if isfield(asi_regimes, name{1})
%             lon = asi_regimes.(name{1});
%             break
%         end
%     end
% 
%     for name = {'lat','latitude','lats'}
%         if isfield(asi_regimes, name{1})
%             lat = asi_regimes.(name{1});
%             break
%         end
%     end
% 
%     lon = lon(:);
%     lat = lat(:);
% 
%     lon0_use = lon0;
%     if max(lon) > 180 && lon0 < 0
%         lon0_use = mod(lon0, 360);
%     elseif max(lon) <= 180 && lon0 > 180
%         lon0_use = lon0 - 360;
%     end
% 
%     P = orient_asi_to_lon_lat_time(param, lon, lat);
% 
%     [~, ilon0] = min(abs(lon - lon0_use));
%     [~, ilat0] = min(abs(lat - lat0));
% 
%     ilon_box = find(abs(lon - lon0_use) <= half_width_deg);
%     ilat_box = find(abs(lat - lat0)     <= half_width_deg);
% 
%     doy = noleap_doy(dates);
% 
%     asi_center = nan(size(doy));
%     asi_std    = nan(size(doy));
% 
%     for ii = 1:numel(doy)
%         if isnan(doy(ii)) || doy(ii) < 1 || doy(ii) > size(P,3)
%             continue
%         end
% 
%         slice = P(:,:,doy(ii));
% 
%         asi_center(ii) = slice(ilon0, ilat0);
% 
%         boxvals = slice(ilon_box, ilat_box);
%         asi_std(ii) = std(boxvals(:), 'omitnan');
%     end
% end

function P = orient_asi_to_lon_lat_time(param, lon, lat)

    sz = size(param);

    if sz(1) == numel(lon) && sz(2) == numel(lat)
        P = param;

    elseif sz(1) == numel(lat) && sz(2) == numel(lon)
        P = permute(param, [2 1 3]);

    elseif sz(2) == numel(lon) && sz(3) == numel(lat)
        P = permute(param, [2 3 1]);

    elseif sz(2) == numel(lat) && sz(3) == numel(lon)
        P = permute(param, [3 2 1]);

    else
        error('Could not infer ASI dimension order. Size is %s', mat2str(sz))
    end
end

function ygrid = grid_by_doy(y, doy, keep, nt)

    ygrid = nan(nt,1);

    y = y(:);
    doy = doy(:);
    keep = keep(:);

    good = keep & isfinite(y) & isfinite(doy) & doy >= 1 & doy <= nt;

    ygrid(doy(good)) = y(good);
end

function shade_between(x, ylo, yhi, color, alpha_val)

    x = x(:);
    ylo = ylo(:);
    yhi = yhi(:);

    good = isfinite(x) & isfinite(ylo) & isfinite(yhi);

    if ~any(good)
        return
    end

    % Break shading across gaps.
    idx = find(good);
    breaks = [1; find(diff(idx) > 1) + 1; numel(idx)+1];

    for bb = 1:numel(breaks)-1
        use_idx = idx(breaks(bb):breaks(bb+1)-1);

        if numel(use_idx) < 2
            continue
        end

        xx = x(use_idx);
        lo = ylo(use_idx);
        hi = yhi(use_idx);

        patch([xx; flipud(xx)], [lo; flipud(hi)], color, ...
              'FaceAlpha', alpha_val, ...
              'EdgeColor', 'none', ...
              'HandleVisibility', 'off');
    end
end

function [lo, hi] = mld_envelope(central, a, b)
% Lower/upper envelope across the central estimate and the two
% MLD-perturbed estimates (Hlo, Hhi), ignoring NaNs.
    lo = min([central, a, b], [], 2, 'omitnan');
    hi = max([central, a, b], [], 2, 'omitnan');
end