% ASI regimes: Example split of one month of data
% Generates Fig. 2 for Kinsella and Mahadevan: "Oceanic processes control air-sea interaction during the1
% South Asian Summer Monsoon"
% Alex Kinsella, February 2026

x = load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/sst_thf_anomaly_corrs_260209.mat');
x = x.sst_thf_anomaly_corrs;
load /Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/oaflux_NIO.mat

%%

lon = 89.5;
lat = 16.5;
dy = 166;

ii = find(flux_turb_NIO.lon==lon,1);
jj = find(flux_turb_NIO.lat==lat,1);

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

% Make plots to show fit
lags = -30:30;
fig1 = makefig;

subplot(2,2,1)
hold on
h1 = plot(lags,SST,'linewidth',3);
h2 = plot(lags,tend,'linewidth',3);
ylim([-0.5 0.5])
ylabel('Correlation','interpreter','latex')
yline(0,'color','k','linewidth',2,'linestyle','--')
xline(0,'color','k','linewidth',2,'linestyle','--')
hold off
makespruce(40)
maketitle('Original Curves (100\%)',40)

subplot(2,2,2)
hold on
plot(lags,c(1)*SST_even_norm,'linewidth',3)
ylim([-0.5 0.5])
plot(lags,c(1)*tend_odd_norm/SST_to_tend,'linewidth',3)
yline(0,'color','k','linewidth',2,'linestyle','--')
xline(0,'color','k','linewidth',2,'linestyle','--')
hold off
makespruce(40)
maketitle(['Ocean-Controlled Curves ($P_o = $',sprintf('%.1f', pct_m1),'\%)'],40)

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
makespruce(40)
maketitle(['Atmosphere-Controlled Curves ($P_a =$',sprintf('%.1f', pct_m2),'\%)'],40)

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
makespruce(40)
maketitle(['Residual (',sprintf('%.1f', pct_res),'\%)'],40)

sgtitle('Ocean-Atmosphere Control Decomposition: Central Bay of Bengal on June 15th','fontsize',50,'interpreter','latex')

disp(['Percent Ocean-Controlled: ',sprintf('%.2f', pct_m1),'%'])
disp(['Percent Atmosphere-Controlled: ',sprintf('%.2f', pct_m2),'%'])
disp(['Percent Residual: ',sprintf('%.2f', pct_res),'%'])
