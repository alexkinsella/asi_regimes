% Perform runs with FCL98 stochastic model
% Save model output to a .mat file
% Alex Kinsella, February 2026

nyears = 100; % Number of years to integrate
ndays = 365 * nyears; % Number of days to integrate
time_days = linspace(0, ndays, ndays); % Time vector for integration

rho0 = 1025; % kg/m^3
cp = 3900; % J/kg/K
%H = 50; % m % FCL98 default
H = 20; % m
%lambda_s = 23.4; % W/m^2/K, FCL98 default
lambda_s = 30; % W/m^2/K
%lambda_r = 1.3; % W/m^2/K, FCL98 default
lambda_r = 5; % W/m^2/K

nu = 1/(rho0*cp*H); % m^2*K/J
alpha = lambda_s*nu; % 1/s
beta = lambda_r*nu; % 1/s
omega_o = 100; % W/m^2

omega_a = linspace(0,10,60);
%omega_a = [0 0.15];

% Noise
tau = 0.1;  % desired timescale (days)
phi = exp(-1/tau);  % AR(1) coefficient
sigma = sqrt(1 - phi^2);

% Initialize array
noise_ocn_AR = zeros(1, ndays);
noise_ocn_AR(1) = 0;  % or another initial condition

% Generate AR(1) in discrete time
for t = 2:ndays
    noise_ocn_AR(t) = phi * noise_ocn_AR(t-1) + sigma * randn;
end

noise_atm = randn(1, ndays); % Stochastic noise for atmosphere
noise_ocn = noise_ocn_AR;
%noise_ocn = randn(1, ndays); % Stochastic noise for ocean

% Interpolate noise for continuous time
noise_atm_interp = @(time) interp1(time_days, noise_atm, time, 'linear', 0);
noise_ocn_interp = @(time) interp1(time_days, noise_ocn, time, 'linear', 0);

% Define the model function
function dTemp_dt = model(time, temps, alpha, beta, nu, omega_a, omega_o,...
                          noise_atm_interp, noise_ocn_interp)
    % temps: ocean temperature anomaly
    SECONDS_PER_DAY = 86400;
    dTemp_dt = -SECONDS_PER_DAY*alpha * (temps - omega_a*noise_atm_interp(time)) - SECONDS_PER_DAY*beta * temps + nu* SECONDS_PER_DAY*omega_o * noise_ocn_interp(time);
end
%%
temp_ocn = nan(numel(time_days)+1,numel(omega_a));

for nn = 1:numel(omega_a)
    % Solve the ODE system using a function handle with additional parameters
    [time_sol, temp_ocn(:,nn)] = ode45(@(time, temps) model(time, temps, alpha, beta, nu, ...
        omega_a(nn), omega_o, noise_atm_interp, noise_ocn_interp), ...
        [0:1:ndays], 0);
    disp(['Done with run ',num2str(nn),' of ',num2str(numel(omega_a))]);
end

Q = lambda_s*(temp_ocn - omega_a.*noise_atm_interp(time_sol));

%% Save
params.rho0 = rho0; 
param.cp = cp;
params.H = H; 
params.lambda_s = lambda_s;
params.lambda_r = lambda_r;
params.omega_o = omega_o;
params.omega_a = omega_a;
params.tau = tau;

FCL.omega_a = omega_a; 
FCL.omega_o = omega_o;
FCL.temp_ocn = temp_ocn;
FCL.Q = Q;
FCL.params = params;
FCL.time_sol = time_sol;
FCL.noise_atm = noise_atm_interp(time_sol);
FCL.noise_ocn = noise_ocn_interp(time_sol);

%save('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/FCL_runs_250114_tropical_params.mat','FCL','-v7.3')
