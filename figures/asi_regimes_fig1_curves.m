% Figure 1: FCL98 stochastic model lagged correlations
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs: P.FCL_runs (pipelines/FCL98_model.m)
% Output: P.fig_dir/asi_regimes_fig1_curves.png (panel plots only; the
%         published Fig. 1 adds annotations in a graphics editor)

clear; clc

P = asi_regimes_paths;

load(P.FCL_runs, 'FCL')

%% Calculate correlation curves
[lagged_corr_temp, lagged_corr_dTemp, lags] = fcl_lagged_corrs(FCL);

%% Plot a single frame
omega_a = FCL.omega_a;

fig1 = makefig(1800, 800);
nn=1;
subplot(1,2,1)
plot(lags, lagged_corr_temp(:,nn), 'LineWidth', 8);
hold on;
plot(lags, lagged_corr_dTemp(:,nn), 'LineWidth', 8,'color','#77AC30');
xlabel('Lag (days)','interpreter','latex');
ylabel('Correlation Coefficient','interpreter','latex');
maketitle({'(a) Ocean Control', '$Q_\mathrm{turb} \sim T_s$'},40);
yline(0,'--','linewidth',2,'color','k','DisplayName','')
grid on;
hold off;
makespruce(40)
ylim([-1 1])
xlim([-30 30])

nn = numel(omega_a);
subplot(1,2,2)
plot(lags, lagged_corr_temp(:,nn), 'LineWidth', 8);
hold on;
plot(lags, lagged_corr_dTemp(:,nn), 'LineWidth', 8,'color','#77AC30');
xlabel('Lag (days)','interpreter','latex');
maketitle({'(b) Atmosphere Control', '$\mathrm{d}T_s/\mathrm{d}t \sim -Q_\mathrm{turb}$'},40);
yline(0,'--','linewidth',2,'color','k','DisplayName','')
legend('$C_{T_s,Q_\mathrm{turb}}(l)$','$C_{\mathrm{d}T_s/\mathrm{d}t, Q_\mathrm{turb}}(l)$','','interpreter','latex','Location', 'northeast');
grid on;
hold off;
makespruce(40)
ylim([-1 1])
xlim([-30 30])

sgtitle('Air-Sea Interaction Regimes: Model Lagged Correlations','interpreter','latex','fontsize',60)

export_figure(fig1, P.fig_dir, 'asi_regimes_fig1_curves')
