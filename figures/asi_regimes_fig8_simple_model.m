% Figure 8: ASI regime parameter P_o in the idealized FCL98 model
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs: P.FCL_runs   (pipelines/FCL98_model.m)
%         P.asi_OAFlux (pipelines/asi_regimes_parameter_calculation.m)
% Output: P.fig_dir/asi_regimes_fig8_simple_model.png

clear; clc

P = asi_regimes_paths;

load(P.FCL_runs, 'FCL')
load(P.asi_OAFlux, 'asi_regimes')

%% Lagged correlations of the model runs

[lagged_corr_temp, lagged_corr_dTemp, lags] = fcl_lagged_corrs(FCL);

omega_a = FCL.omega_a;
omega_o = FCL.omega_o;
omega_a_rescale = FCL.params.lambda_s*omega_a;

%% Apply decomposition algorithm (unweighted) to each model run

c = nan([size(lagged_corr_temp,2),2]);

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

    c(nn,:) = lsqnonneg(C,d);
end

param =  c(:,1).^2./(c(:,1).^2+c(:,2).^2);

%% 5th-95th percentile range of the monthly climatological OAFlux P_o

% Calendar month of each day in the 365-day no-leap year
month_of_doy = month(datetime(2021,1,1:365));

param_monthly = nan([size(asi_regimes.param,[1 2]), 12]);
for mm = 1:12
    param_monthly(:,:,mm) = mean(asi_regimes.param(:,:,month_of_doy==mm), 3, 'omitnan');
end

asi_5th = prctile(param_monthly(:), 5);
asi_95th = prctile(param_monthly(:), 95);
fprintf('5th-95th percentile of monthly climatological P_o: %.4f to %.4f\n', asi_5th, asi_95th)

% Corresponding range of forcing ratios, by inverting P_o(ratio)
x_ratio = omega_o./omega_a_rescale;
x_ratio = x_ratio(2:end);

ratio_5th = interp1(param(2:end),x_ratio,asi_5th);
ratio_95th = interp1(param(2:end),x_ratio,asi_95th);

%% Plot param vs hypothesis

fig1 = makefig(1800, 987);  % window size of the published figure
plot(omega_o./omega_a_rescale,param,'linewidth',2)
hold on
plot(omega_o./omega_a_rescale,(omega_o./sqrt(omega_a_rescale.^2+omega_o.^2)).^2,'linewidth',2)
ylabel('Parameter Value','interpreter','latex')
xregion(ratio_5th,ratio_95th)

xlabel('Ocean/Atmosphere Forcing Ratio $\omega_o / ( \lambda_s \omega_a )$','interpreter','latex')
legend('ASI Regime Parameter $P_o$',sprintf('%s\n%s', 'Fraction Ocean Forcing', '$(\omega_o)^2/((\lambda_s\omega_a)^2+\omega_o^2)$'),'NIO 90$\%$ interpercentile range','interpreter','latex','location','southeast')
xlim([0 4])
makespruce(30)
maketitle({'ASI Regime Parameter $P_o$ in an Idealized Stochastic Model'},40)

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

export_figure(fig1, P.fig_dir, 'asi_regimes_fig8_simple_model')
