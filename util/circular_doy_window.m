function inds = circular_doy_window(doy365, center_doy, halfwidth)
%CIRCULAR_DOY_WINDOW Logical index for circular day-of-year window.
%
% Alex Kinsella, June 2026

dist = abs(doy365 - center_doy);
dist = min(dist, 365 - dist);
inds = dist <= halfwidth;
end

