% Create daily climatology from spatially-distributed time series
function varclimday = clim_daily(data,dvec,nclim)
% Data is time series, dvec is associated datevec for measurements, and
% nclim is number of days on either side to average over.
% Time should be final dimension of data array.
% For example, data may be OISST.SST, dvec is OISST.dvec, nclim is 2.
%
% Alex Kinsella, July 28, 2022

sz = size(data);
if sum(sz) == prod(sz) + 1 % Checking for time series at a point
    if sz(1)>sz(2)
        sz = flip(sz); % Putting time coord last
    end
end
datavec = reshape(data,[numel(data)/sz(end) sz(end)]);
varclimdayvec = nan([numel(data)/sz(end) 366]);
for dd = 1:366
    ddinds = ismember(day(datetime(dvec),'dayofyear'),dd-nclim:dd+nclim);
    varclimdayvec(:,dd) = nanmean(datavec(:,ddinds),2);
end

varclimday = reshape(varclimdayvec,[sz(1:end-1), 366]);

end