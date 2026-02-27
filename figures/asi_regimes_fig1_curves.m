% Make fig 1 from FCL98 stochastic model
% Alex Kinsella, February 2026

load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/FCL_runs_250114_tropical_params.mat')

%% Calculate correlation curves with SHF
Q = FCL.Q;
temp_ocn = FCL.temp_ocn;
time_sol = FCL.time_sol;

surface_heat_flux = Q;

% Time step (assuming uniform spacing from ode45)
dt = time_sol(2) - time_sol(1);

% Calculate dT_ocn/dt using central differences
dTemp_ocn_dt = zeros(size(temp_ocn));
dTemp_ocn_dt(2:end-1,:) = (temp_ocn(3:end,:) - temp_ocn(1:end-2,:)) / (2 * dt); % Central difference
dTemp_ocn_dt(1,:) = (temp_ocn(2,:) - temp_ocn(1,:)) / dt; % Forward difference for the first point
dTemp_ocn_dt(end,:) = (temp_ocn(end,:) - temp_ocn(end-1,:)) / dt; % Backward difference for the last point

% Define lag range
max_lag = 30; % Maximum lag in days
lags = -max_lag:max_lag;

% Initialize arrays for lagged correlation coefficients
lagged_corr_temp = zeros(length(lags), size(temp_ocn,2)); % Corr(temp_ocn, surface_heat_flux)
lagged_corr_dTemp = zeros(length(lags), size(temp_ocn,2)); % Corr(dTemp_ocn/dt, surface_heat_flux)

% Compute lagged correlations
for nn = 1:size(temp_ocn,2)
    for i = 1:length(lags)
        lag = lags(i);
        if lag < 0
            % Negative lag: shift series backward
            lagged_corr_temp(i,nn) = corr(temp_ocn(1:end+lag,nn), surface_heat_flux(1-lag:end,nn));
            lagged_corr_dTemp(i,nn) = corr(dTemp_ocn_dt(1:end+lag,nn), surface_heat_flux(1-lag:end,nn));
        elseif lag > 0
            % Positive lag: shift series forward
            lagged_corr_temp(i,nn) = corr(temp_ocn(1+lag:end,nn), surface_heat_flux(1:end-lag,nn));
            lagged_corr_dTemp(i,nn) = corr(dTemp_ocn_dt(1+lag:end,nn), surface_heat_flux(1:end-lag,nn));
        else
            % Zero lag: no shift
            lagged_corr_temp(i,nn) = corr(temp_ocn(:,nn), surface_heat_flux(:,nn));
            lagged_corr_dTemp(i,nn) = corr(dTemp_ocn_dt(:,nn), surface_heat_flux(:,nn));
        end
    end
end

%% Plot a single frame
omega_a = FCL.omega_a;
omega_o = FCL.omega_o;

fig1 = figure('color',[1 1 1],'position',[100 100 1800 800]);
nn=1;
subplot(1,2,1)
plot(lags, lagged_corr_temp(:,nn), 'LineWidth', 8);
hold on;
plot(lags, lagged_corr_dTemp(:,nn), 'LineWidth', 8,'color','#77AC30');
xlabel('Lag (days)','interpreter','latex');
ylabel('Correlation Coefficient','interpreter','latex');
maketitle('Ocean Control: $Q_\mathrm{turb} \sim T_s$',40);
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
maketitle('Atmosphere Control: $\mathrm{d}T_s/\mathrm{d}t \sim -Q_\mathrm{turb}$',40);
yline(0,'--','linewidth',2,'color','k','DisplayName','')
legend('$C_{T_s,Q_\mathrm{turb}}(l)$','$C_{\mathrm{d}T_s/\mathrm{d}t, Q_\mathrm{turb}}(l)$','','interpreter','latex','Location', 'northeast');
grid on;
hold off;
makespruce(40)
ylim([-1 1])
xlim([-30 30])

sgtitle('Air-Sea Interaction Regimes: Model Lagged Correlations','interpreter','latex','fontsize',60)