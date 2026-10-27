% Process daily CERES SYN1deg all-sky surface radiative fluxes over the
% northern Indian Ocean (Section 4b).
%
% Kinsella, A., and A. Mahadevan, 2026: A New Diagnostic of Air-Sea
% Interaction Reveals Ocean Control of Turbulent Heat Flux in the South
% Asian Summer Monsoon. J. Climate, 39, 6319-6333,
% https://doi.org/10.1175/JCLI-D-26-0118.1.
%
% Inputs:  P.ceres_file, CERES SYN1deg-Day Ed4.1 subset (40-100E,
%          10S-30N, March 2000 - March 2024) with adjusted all-sky
%          surface SW and LW up/down fluxes
% Outputs: P.ceres_rad, structure CERES_rad with fields lon, lat, dnum
%          (datenum), dvec, SWnet and LWnet (W/m^2, positive downward)

clear; clc

P = asi_regimes_paths;

ncpath_flux = P.ceres_file;

LWup = ncread(ncpath_flux,'adj_atmos_lw_up_all_surface_daily');
LWdwn = ncread(ncpath_flux,'adj_atmos_lw_down_all_surface_daily');
SWup = ncread(ncpath_flux,'adj_atmos_sw_up_all_surface_daily');
SWdwn = ncread(ncpath_flux,'adj_atmos_sw_down_all_surface_daily');

lon = ncread(ncpath_flux,'lon');
lat = ncread(ncpath_flux,'lat');
dnum = datenum(datetime(2000,3,1)+days(ncread(ncpath_flux,'time'))); % time: days since 2000-03-01
dvec = datevec(dnum);

CERES_rad.dnum = dnum;
CERES_rad.dvec = dvec;
CERES_rad.lon = lon;
CERES_rad.lat = lat;
CERES_rad.LWnet = LWdwn-LWup;
CERES_rad.SWnet = SWdwn-SWup;

if ~exist(P.output_dir,'dir')
    mkdir(P.output_dir)
end
save(P.ceres_rad,'CERES_rad','-v7.3')
