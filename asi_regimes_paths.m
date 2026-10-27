function P = asi_regimes_paths()
% Paths to input datasets and intermediate files for the asi_regimes code.
% Edit the USER SETTINGS section to match your system; every script in
% this repository reads its file locations from here.
%
% Usage (from any script):
%   P = asi_regimes_paths;
%
% Calling this function also adds the repository's util/ and util/gsw/
% folders to the MATLAB path.
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.

%% ---------------- USER SETTINGS ----------------

% Input datasets (see README for download sources)
P.oaflux_dir  = '/Volumes/proj/mahadevanlab/datasets/OAFlux/daily/';                % OAFlux daily NetCDF, with turbulence/, radiation/, ... subfolders
P.era5_dir    = '/Volumes/proj/mahadevanlab/datasets/ERA5/Indian_Ocean/THF_SST/';   % ERA5 THF + SST NetCDF files, one per year
P.data_dir    = '/Volumes/proj/mahadevanlab/kinsella/data/asi_regimes/';           % CERES subset and mooring files

% Where intermediate .mat files are written and read
P.output_dir  = '/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes_release/';

% Where figures and the supplemental video are written
P.fig_dir     = fullfile(fileparts(mfilename('fullpath')), 'output_figures');

%% ---------------- DERIVED PATHS (no need to edit) ----------------

% CERES SYN1deg daily surface fluxes
P.ceres_file = fullfile(P.data_dir, 'CERES_SYN1deg-Day_Terra-Aqua-MODIS_Ed4.1_Subset_20000301-20240331.nc');

% Arabian Sea Experiment mooring (EPIC files)
P.AS_flux  = fullfile(P.data_dir, 'arabianseaflux_1hr.epic');
P.AS_temp  = fullfile(P.data_dir, 'arabianseatemp.nc.epic');
P.AS_sal   = fullfile(P.data_dir, 'arabianseasal.nc.epic');
P.AS_curr  = fullfile(P.data_dir, 'arabianseavel.epic');

% ASIRI Bay of Bengal mooring
P.BoB_flux = fullfile(P.data_dir, 'beng1_flux_1hr.mat');
P.BoB_ts   = fullfile(P.data_dir, 'bob1_zsal_hr_vfilled.mat');
P.BoB_curr = fullfile(P.data_dir, 'bob1_zcur_hr.mat');

% Natural Earth coastlines (public domain, shipped in util/)
util_dir = fullfile(fileparts(mfilename('fullpath')), 'util');
P.coast_10m = fullfile(util_dir, 'ne_10m_coastline', 'ne_10m_coastline.shp');
P.coast_50m = fullfile(util_dir, 'ne_50m_coastline', 'ne_50m_coastline.shp');

% Intermediate files, in the order the pipelines produce them
P.oaflux_NIO       = fullfile(P.output_dir, 'oaflux_NIO.mat');                % pipelines/asi_regimes_oaflux_processing.m
P.era5_1985_2020   = fullfile(P.output_dir, 'ERA5_THF_SST_IO_1985_2020.mat'); % pipelines/asi_regimes_era5_processing.m
P.FCL_runs         = fullfile(P.output_dir, 'FCL_runs_tropical_params.mat');  % pipelines/FCL98_model.m
P.corrs_OAFlux     = fullfile(P.output_dir, 'sst_thf_anomaly_corrs_OAFlux.mat');  % pipelines/asi_regimes_parameter_calculation.m
P.asi_OAFlux       = fullfile(P.output_dir, 'asi_regimes_OAFlux.mat');            % pipelines/asi_regimes_parameter_calculation.m
P.corrs_ERA5       = fullfile(P.output_dir, 'sst_thf_anomaly_corrs_ERA5.mat');    % pipelines/asi_regimes_parameter_calculation_ERA5.m
P.asi_ERA5         = fullfile(P.output_dir, 'asi_regimes_ERA5.mat');              % pipelines/asi_regimes_parameter_calculation_ERA5.m
P.monte_carlo      = fullfile(P.output_dir, 'sst_thf_anomaly_corrs_monte_carlo.mat'); % pipelines/monte_carlo_resampling.m
P.ceres_rad        = fullfile(P.output_dir, 'CERES_rad.mat');                 % pipelines/asi_regimes_ceres_processing.m

%% ---------------- SETUP ----------------

addpath(util_dir)
addpath(fullfile(util_dir, 'gsw'))  % subset of the GSW (TEOS-10) toolbox, v3.06.13

if ~exist(P.output_dir, 'dir')
    warning('asi_regimes_paths:noOutputDir', ...
        'Output directory does not exist yet: %s', P.output_dir)
end

end
