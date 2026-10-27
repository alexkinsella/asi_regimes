function [doy365, keep, dt_keep] = doy_noleap(dt)
%DOY_NOLEAP Return 365-day no-leap day-of-year index.
%
% Feb. 29 is omitted. Dates after Feb. 28 in leap years are shifted back
% by one day, so March 1 is always day 60 and Dec. 31 is always day 365.
%
% Alex Kinsella, June 2026

dt = dt(:);
keep = ~(month(dt)==2 & day(dt)==29);
dt_keep = dt(keep);

doy = day(dt_keep,'dayofyear');
is_leap_year = eomday(year(dt_keep),2)==29;
after_feb28 = is_leap_year & month(dt_keep) > 2;

doy365 = doy - after_feb28;
end

