% Figure 2: Example split into ASI regimes
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Unweighted decomposition of the lagged correlation curves at one grid
% point and day; panel titles give each part's share of the variance.
%
% Inputs: P.corrs_OAFlux, P.asi_OAFlux
%         (pipelines/asi_regimes_parameter_calculation.m)
% Output: P.fig_dir/asi_regimes_fig2_decomposition.png

clear; clc

P = asi_regimes_paths;

x = load(P.corrs_OAFlux, 'sst_thf_anomaly_corrs');
x = x.sst_thf_anomaly_corrs;
load(P.asi_OAFlux, 'asi_regimes')   % for the lon/lat grid

%% Process data

lon = 89.5;
lat = 16.5;
dy = 166;   % June 15 in the 365-day no-leap calendar

ii = find(asi_regimes.lon==lon,1);
jj = find(asi_regimes.lat==lat,1);

SST = squeeze(x.corrs_SST(ii,jj,:,dy));
tend = squeeze(x.corrs_tend(ii,jj,:,dy));

SST_even = 1/2*(SST+flip(SST));
SST_odd = 1/2*(SST-flip(SST));
tend_even = 1/2*(tend+flip(tend));
tend_odd = 1/2*(tend-flip(tend));

SST_even_norm = SST_even/norm(SST_even);
SST_odd_norm = SST_odd/norm(SST_odd);
tend_even_norm = tend_even/norm(tend_even);
tend_odd_norm = tend_odd/norm(tend_odd);

SST_to_tend = norm(SST)/norm(tend);

C = [[SST_even_norm; tend_odd_norm./SST_to_tend], [SST_odd_norm; tend_even_norm./SST_to_tend]];
d = [SST;tend];

[c,resnorm,residual] = lsqnonneg(C,d);

% --- Data-space decomposition ---
m1 = c(1) * C(:,1);
m2 = c(2) * C(:,2);
model = m1 + m2;

% Variance fractions 
pct_m1   = 100 * dot(m1, model) / norm(d)^2;
pct_m2   = 100 * dot(m2, model) / norm(d)^2;
pct_res  = 100 * norm(residual)^2 / norm(d)^2;

%% Make plots to show fit
lags = -30:30;
fig1 = makefig(1800, 987);  % window size of the published figure

subplot(2,2,1)
hold on
h1 = plot(lags,SST,'linewidth',3);
h2 = plot(lags,tend,'linewidth',3);
ylim([-0.5 0.5])
ylabel('Correlation','interpreter','latex')
yline(0,'color','k','linewidth',2,'linestyle','--')
xline(0,'color','k','linewidth',2,'linestyle','--')
hold off
makespruce(30)
maketitle('(a) Original Curves (100\%)',30)

subplot(2,2,2)
hold on
plot(lags,c(1)*SST_even_norm,'linewidth',3)
ylim([-0.5 0.5])
plot(lags,c(1)*tend_odd_norm/SST_to_tend,'linewidth',3)
yline(0,'color','k','linewidth',2,'linestyle','--')
xline(0,'color','k','linewidth',2,'linestyle','--')
hold off
makespruce(30)
maketitle(['(b) Ocean-Controlled Curves ($P_o = $',sprintf('%.1f', pct_m1),'\%)'],30)

subplot(2,2,3)
hold on
plot(lags,c(2)*SST_odd_norm,'linewidth',3)
ylim([-0.5 0.5])
plot(lags,c(2)*tend_even_norm/SST_to_tend,'linewidth',3)
yline(0,'color','k','linewidth',2,'linestyle','--')
xline(0,'color','k','linewidth',2,'linestyle','--')
hold off
ylabel('Correlation','interpreter','latex')
xlabel('Lag (Days)','interpreter','latex')
makespruce(30)
maketitle(['(c) Atmosphere-Controlled Curves ($P_a =$',sprintf('%.1f', pct_m2),'\%)'],30)

subplot(2,2,4)
hold on
plot(lags,residual(1:61),'linewidth',3)
ylim([-0.5 0.5])
plot(lags,residual(62:end),'linewidth',3)
yline(0,'color','k','linewidth',2,'linestyle','--')
xline(0,'color','k','linewidth',2,'linestyle','--')
hold off
xlabel('Lag (Days)','interpreter','latex')
legend('$T_s$-$Q_\mathrm{turb}$ Corr','$\mathrm{d}T_s/\mathrm{d}t$-$Q_\mathrm{turb}$ Corr','interpreter','latex','location','northeast')
makespruce(30)
maketitle(['(d) Residual (',sprintf('%.1f', pct_res),'\%)'],30)

sgtitle('Ocean-Atmosphere Control Decomposition: Central Bay of Bengal on June 15th','fontsize',40,'interpreter','latex')

disp(['Percent Ocean-Controlled: ',sprintf('%.2f', pct_m1),'%'])
disp(['Percent Atmosphere-Controlled: ',sprintf('%.2f', pct_m2),'%'])
disp(['Percent Residual: ',sprintf('%.2f', pct_res),'%'])

export_figure(fig1, P.fig_dir, 'asi_regimes_fig2_decomposition')
