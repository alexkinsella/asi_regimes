function out = circular_movmean_doy(data, win)
%CIRCULAR_MOVMEAN_DOY Circular moving mean along the day-of-year dimension.
%
% data is a 3-D array (lon x lat x 365) with day of year as the third
% dimension. win is the window length in days and must be odd.
%
% Alex Kinsella, June 2026

n = size(data, ndims(data));

if mod(win,2)==0
    error('Use an odd smoothing window for circular_movmean_doy.')
end

halfwin = floor(win/2);

data_ext = cat(ndims(data), ...
               data(:,:,n-halfwin+1:n), ...
               data, ...
               data(:,:,1:halfwin));

out_ext = movmean(data_ext, win, ndims(data_ext), 'omitnan');
out = out_ext(:,:,halfwin+1:halfwin+n);
end
