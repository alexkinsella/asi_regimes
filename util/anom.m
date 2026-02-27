% Create anom from clim_daily

function [data_anom,data_clim] = anom(data,dvec,nclim)

sz = size(data);
if sum(sz) == prod(sz) + 1 % Checking for time series at a point
    if sz(1)>sz(2)
        data = data'; % Putting time coord last
    end
end

data_clim = clim_daily(data,dvec,nclim);
data_anom = nan(size(data));

if numel(size(data))==3
    for tt = 1:size(dvec,1)
        dd = day(datetime(dvec(tt,:)),'dayofyear');
        data_anom(:,:,tt) = data(:,:,tt) - data_clim(:,:,dd);
    end
elseif numel(size(data))==2
        for tt = 1:size(dvec,1)
        dd = day(datetime(dvec(tt,:)),'dayofyear');
        data_anom(:,tt) = data(:,tt) - data_clim(:,dd);
        end
end

end