function varclimday = clim_daily_noleap(data,doy365,nclim)
%CLIM_DAILY_NOLEAP Daily climatology on circular 365-day calendar.
%
% data has time as final dimension.
% doy365 is a 365-day no-leap day-of-year index.
% nclim is number of days on either side to average over.
%
% Alex Kinsella, June 2026

sz = size(data);

if sum(sz) == prod(sz) + 1
    if sz(1)>sz(2)
        sz = flip(sz);
    end
end

datavec = reshape(data,[numel(data)/sz(end) sz(end)]);
varclimdayvec = nan([numel(data)/sz(end) 365]);

for dd = 1:365
    ddinds = circular_doy_window(doy365, dd, nclim);
    varclimdayvec(:,dd) = mean(datavec(:,ddinds),2,'omitnan');
end

varclimday = reshape(varclimdayvec,[sz(1:end-1), 365]);
end

