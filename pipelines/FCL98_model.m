% Run the FCL98 stochastic mixed layer model (Section 2a, Eqs. 2-3,
% Table 1) for a range of atmospheric forcing strengths. Used in Figs. 1
% and 8.
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Outputs: P.FCL_runs, structure FCL with the ocean temperature anomaly
%          temp_ocn (K) and turbulent heat flux Q (W/m^2, positive
%          upward) for each atmospheric forcing strength omega_a (K)
%
% The model is driven by random noise. rng(seed) makes a rerun
% reproducible; the published figures used an earlier, unseeded run.

clear; clc

P = asi_regimes_paths;

seed = 0;
rng(seed, 'twister')

%% Parameters (Table 1)

nyears = 100; % Number of years to integrate
ndays = 365 * nyears; % Number of days to integrate
time_days = linspace(0, ndays, ndays); % Time vector for noise

rho0 = 1025; % kg/m^3
cp = 3900; % J/kg/K
H = 20; % m (FCL98: 50 m)
lambda_s = 30; % W/m^2/K, turbulent feedback (FCL98: 23.4)
lambda_r = 5; % W/m^2/K, radiative feedback (FCL98: 1.3)

nu = 1/(rho0*cp*H); % m^2*K/J
alpha = lambda_s*nu; % 1/s
beta = lambda_r*nu; % 1/s
omega_o = 100; % W/m^2, oceanic forcing strength

omega_a = linspace(0,10,60); % K, atmospheric forcing strengths

%% Stochastic forcing

% Oceanic forcing: AR(1) process with memory time scale tau (Eq. 3)
tau = 0.1;  % days
phi = exp(-1/tau);
sigma = sqrt(1 - phi^2);

noise_ocn = zeros(1, ndays);
for t = 2:ndays
    noise_ocn(t) = phi * noise_ocn(t-1) + sigma * randn;
end

% Atmospheric forcing: white noise
noise_atm = randn(1, ndays);

% Interpolate noise for continuous time
noise_atm_interp = @(time) interp1(time_days, noise_atm, time, 'linear', 0);
noise_ocn_interp = @(time) interp1(time_days, noise_ocn, time, 'linear', 0);

%% Integrate the model for each atmospheric forcing strength

temp_ocn = nan(numel(time_days)+1,numel(omega_a));

for nn = 1:numel(omega_a)
    [time_sol, temp_ocn(:,nn)] = ode45(@(time, temps) model(time, temps, alpha, beta, nu, ...
        omega_a(nn), omega_o, noise_atm_interp, noise_ocn_interp), ...
        0:1:ndays, 0);
    fprintf('Done with run %d of %d\n', nn, numel(omega_a));
end

% Turbulent heat flux, positive upward
Q = lambda_s*(temp_ocn - omega_a.*noise_atm_interp(time_sol));

%% Save

params.rho0 = rho0;
params.cp = cp;
params.H = H;
params.lambda_s = lambda_s;
params.lambda_r = lambda_r;
params.omega_o = omega_o;
params.omega_a = omega_a;
params.tau = tau;
params.seed = seed;

FCL.omega_a = omega_a;
FCL.omega_o = omega_o;
FCL.temp_ocn = temp_ocn;
FCL.Q = Q;
FCL.params = params;
FCL.time_sol = time_sol;
FCL.noise_atm = noise_atm_interp(time_sol);
FCL.noise_ocn = noise_ocn_interp(time_sol);

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.FCL_runs, 'FCL', '-v7.3')

%% Model (Eq. 2)

function dTemp_dt = model(time, temps, alpha, beta, nu, omega_a, omega_o, ...
                          noise_atm_interp, noise_ocn_interp)
    % temps: ocean temperature anomaly (K); time in days
    SECONDS_PER_DAY = 86400;
    dTemp_dt = -SECONDS_PER_DAY*alpha * (temps - omega_a*noise_atm_interp(time)) ...
               - SECONDS_PER_DAY*beta * temps ...
               + nu*SECONDS_PER_DAY*omega_o * noise_ocn_interp(time);
end
