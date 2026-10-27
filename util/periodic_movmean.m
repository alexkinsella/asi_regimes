function movmean_result = periodic_movmean(data, window)
    % Create a periodic extension of data for moving mean calculation
    extended_data = [data(end - window + 1:end); data; data(1:window)];
    movmean_extended = movmean(extended_data, window);
    movmean_result = movmean_extended(window + 1:end - window);
end