% Figure to show use of ASI regime parameter in models
% Alex Kinsella, February 2026

load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/FCL_runs_250114_tropical_params.mat')
load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_OAFlux_260210.mat')

%%
Q = FCL.Q;
temp_ocn = FCL.temp_ocn;
time_sol = FCL.time_sol;

surface_heat_flux = Q;

% Time step 
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

%% Find threshold for correlation amplitudes
omega_a = FCL.omega_a;
omega_o = FCL.omega_o;
omega_a_rescale = FCL.params.lambda_s*omega_a;

for nn = 1:numel(omega_a)
    corr_mag(nn) = max(abs(lagged_corr_temp(:,nn)))+max(abs(lagged_corr_dTemp(:,nn)));
end

%% Apply decomposition algorithm 

lags = -30:30;
c = nan([size(lagged_corr_temp,2),2]);
resnorm = nan(size(lagged_corr_temp,2),1);
residual = nan(size(lagged_corr_temp,2),2*size(lagged_corr_temp,1));

for nn = 1:size(lagged_corr_temp,2)
    SST = lagged_corr_temp(:,nn);
    tend = lagged_corr_dTemp(:,nn);
    
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
    
    [c(nn,:),resnorm(nn),residual(nn,:)] = lsqnonneg(C,d);
end

param =  c(:,1).^2./(c(:,1).^2+c(:,2).^2);


%% Calculate range for NIO
x_ratio = omega_o./omega_a_rescale;
x_ratio = x_ratio(2:end);

asi_95th = 0.9486; % 95th percentile
asi_5th = 0.1380; % 5th percentile

ratio_5th = interp1(param(2:end),x_ratio,asi_5th);
ratio_95th = interp1(param(2:end),x_ratio,asi_95th);

%% Make histogram datasets
% Define common bins
edges = 0:0.05:1;

% Datasets
param_summer = asi_regimes.param(:,:,152:243); % JJA
param_spring = asi_regimes.param(:,:,60:151); % MAM

% Histogram counts
c1 = histcounts(param_summer, edges);
c2 = histcounts(param_spring, edges);

c1 = c1 / sum(c1);
c2 = c2 / sum(c2);

% Bin centers
centers = edges(1:end-1) + diff(edges)/2;
bw = diff(edges(1:2));   % bin width

% Offset for side-by-side bars
offset = 0.35 * bw;

%% Plot param vs hypothesis
omega_a = FCL.omega_a;
omega_o = FCL.omega_o;
omega_a_rescale = FCL.params.lambda_s*FCL.omega_a;

fig1 = makefig;
plot(omega_o./omega_a_rescale,param,'linewidth',2)
hold on
plot(omega_o./omega_a_rescale,(omega_o./sqrt(omega_a_rescale.^2+omega_o.^2)).^2,'linewidth',2)
ylabel('Parameter Value','interpreter','latex')
xregion(ratio_5th,ratio_95th)

xlabel('Ocean/Atmosphere Forcing Ratio $\omega_o / ( \lambda_s \omega_a )$','interpreter','latex')
legend('ASI Regime Parameter $P_o$',sprintf('%s\n%s', 'Percentage Ocean Forcing', '$(\omega_o)^2/((\lambda_s\omega_a)^2+\omega_o^2)$'),'NIO 90$\%$ interpercentile range','interpreter','latex','location','southeast')
xlim([0 4])
makespruce(30)
maketitle({'ASI Regime Parameter $P_o$ in a Simple Model'},40)

% Near y = 1 (top)
text( ...
    -0.09, 0.93, sprintf('%s\n%s','Ocean','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'FontSize',25, ...
    'Interpreter','latex','rotation', 90)

% Near y = 0 (bottom)
text( ...
    -0.04, -0.03, sprintf('%s\n%s','Atm.','control'), ...
    'Units','normalized', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','bottom', ...
    'FontSize',25, ...
    'Interpreter','latex','Rotation', 90)

