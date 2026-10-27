function [lagged_corr_temp, lagged_corr_dTemp, lags] = fcl_lagged_corrs(FCL)
%FCL_LAGGED_CORRS Lagged correlations from the FCL98 stochastic model runs.
%
% [lagged_corr_temp, lagged_corr_dTemp, lags] = fcl_lagged_corrs(FCL)
%
% For each model run (column of FCL.temp_ocn), correlates the ocean
% temperature anomaly and its tendency with the turbulent heat flux FCL.Q
% at lags of -30 to 30 days (positive lag: temperature lags Q).
%
% Alex Kinsella, October 2026

Q = FCL.Q;
temp_ocn = FCL.temp_ocn;
time_sol = FCL.time_sol;

surface_heat_flux = Q;

% Time step (uniform output times from ode45)
dt = time_sol(2) - time_sol(1);

% Calculate dT_ocn/dt using central differences
dTemp_ocn_dt = zeros(size(temp_ocn));
dTemp_ocn_dt(2:end-1,:) = (temp_ocn(3:end,:) - temp_ocn(1:end-2,:)) / (2 * dt); % Central difference
dTemp_ocn_dt(1,:) = (temp_ocn(2,:) - temp_ocn(1,:)) / dt; % Forward difference for the first point
dTemp_ocn_dt(end,:) = (temp_ocn(end,:) - temp_ocn(end-1,:)) / dt; % Backward difference for the last point

% Define lag range
max_lag = 30; % Maximum lag in days
lags = -max_lag:max_lag;

lagged_corr_temp = zeros(length(lags), size(temp_ocn,2)); % Corr(temp_ocn, surface_heat_flux)
lagged_corr_dTemp = zeros(length(lags), size(temp_ocn,2)); % Corr(dTemp_ocn/dt, surface_heat_flux)

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

end
