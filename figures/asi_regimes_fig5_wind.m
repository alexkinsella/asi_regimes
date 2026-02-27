% Makes wind figure for ASI regimes manuscript
% Alex Kinsella, February 2026

%% Load data
load('/Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/asi_regimes_weighted_OAFlux_260210.mat');
load /Volumes/proj/mahadevanlab/kinsella/output/asi_regimes/manuscript/oaflux_NIO.mat

%% Make OAFlux wind clim
nhalf = 5;
nsmooth = 30;

% Compute OAflux wind clim
ws_clim_OAflux = clim_daily(flux_turb_NIO.ws, datevec(flux_turb_NIO.time), nhalf);
ws_clim_OAflux = movmean(ws_clim_OAflux, nsmooth, 3);
ws_clim_OAflux_short = ws_clim_OAflux(:,:,1:365);

%% Combine Arabian Sea and Bay of Bengal

lonmatrix = repmat(asi_regimes.lon,[1 size(asi_regimes.param,[2 3])]);
latmatrix = permute(repmat(asi_regimes.lat,[1 size(asi_regimes.param,[1 3])]),[2 1 3]);

% Line from (5N, 50E) to (26N, 62E)
m = (62 - 50) / (26 - 5);  % slope
b = 50 - m * 5;            % intercept

% Arabian Sea
reg1_lonmin = 50;
reg1_lonmax = 80;
reg1_latmin = 5;
reg1_latmax = 26;

% Bay of Bengal
reg2_lonmin = 80;
reg2_lonmax = 98;
reg2_latmin = 5;
reg2_latmax = 23;

reg1_mask = ...
    lonmatrix >= 50 & lonmatrix <= 80 & ...
    latmatrix >= 5  & latmatrix <= 26 & ...
    lonmatrix >= (m * latmatrix + b);  % Arabian Sea mask

% Define spatial mask for Bay of Bengal
reg2_mask = ...
    lonmatrix >= 80 & lonmatrix <= 98 & ...
    latmatrix >= 5  & latmatrix <= 23;

% Define AS + BoB polygon outlines for plotting

% Arabian Sea polygon (rectangle clipped by the line lon = m*lat + b)
AS_lon = [50 80 80 62 50];
AS_lat = [ 5  5 26 26  5];

% Bay of Bengal rectangle
BoB_lon = [80 98 98 80 80];
BoB_lat = [ 5  5 23 23  5];

% Apply masks to 3D arrays

ws_as = ws_clim_OAflux_short(reg1_mask);
ws_bob = ws_clim_OAflux_short(reg2_mask);
ws_both = ws_clim_OAflux_short(reg1_mask|reg2_mask);

asi_as  = asi_regimes.param(reg1_mask);
asi_bob = asi_regimes.param(reg2_mask);
asi_both = asi_regimes.param(reg1_mask|reg2_mask);


x_bob = ws_bob(:);
x_as = ws_as(:);
x_both = ws_both(:);

y_bob = asi_bob(:);
y_as = asi_as(:);
y_both = asi_both(:);

% Create 2D histogram
numBins = 100;  % Set the number of bins
[counts_bob, edgesX_bob, edgesY_bob] = histcounts2(x_bob, y_bob, numBins);
[counts_as, edgesX_as, edgesY_as] = histcounts2(x_as, y_as, numBins);
[counts_both, edgesX_both, edgesY_both] = histcounts2(x_both, y_both, numBins);

% Define bin centers for plotting
binCentersX_bob = (edgesX_bob(1:end-1) + edgesX_bob(2:end)) / 2;
binCentersY_bob = (edgesY_bob(1:end-1) + edgesY_bob(2:end)) / 2;
binCentersX_as = (edgesX_as(1:end-1) + edgesX_as(2:end)) / 2;
binCentersY_as = (edgesY_as(1:end-1) + edgesY_as(2:end)) / 2;
binCentersX_both = (edgesX_both(1:end-1) + edgesX_both(2:end)) / 2;
binCentersY_both = (edgesY_both(1:end-1) + edgesY_both(2:end)) / 2;

% Calculate mean Y value in each X bin
meanYInXBin_bob = nan(1, numBins);
meanYInXBin_as = nan(1, numBins);
meanYInXBin_both = nan(1, numBins);

for i = 1:numBins
    % Find indices where x falls in the current bin
    inXBin_bob = x_bob >= edgesX_bob(i) & x_bob < edgesX_bob(i+1);
    inXBin_as = x_as >= edgesX_as(i) & x_as < edgesX_as(i+1);
    inXBin_both = x_both >= edgesX_both(i) & x_both < edgesX_both(i+1);

    % Calculate mean Y value for those points
    if any(inXBin_bob)
        meanYInXBin_bob(i) = nanmedian(y_bob(inXBin_bob));
    end
    if any(inXBin_as)
        meanYInXBin_as(i) = nanmedian(y_as(inXBin_as));
    end
    if any(inXBin_both)
        meanYInXBin_both(i) = nanmedian(y_both(inXBin_both));
    end
end

% Plot the density as a 2D image
fig1=makefig;

ax1 = subplot(2,2,3);
imagesc(binCentersX_both, binCentersY_both, (counts_both.'));
axis xy;  
cb = makecb([0 100],cmocean('amp'),'Count',30,ax1);
xlabel('Windspeed','interpreter','latex');
ylabel('ASI Parameter $P_o$','interpreter','latex');
xlim([2.7 12])
makespruce(30)
maketitle('(c) Windspeed-ASI Regime Density Plot',40);

% Overlay mean Y values for each X bin as red dots
hold on;
for i = 1:numBins
    if ~isnan(meanYInXBin_as(i))
plot(binCentersX_both(i), meanYInXBin_both(i), 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'r');    
    end
end
hold off;

% Now make maps

land_mask = double(isnan(asi_regimes.param(:,:,1)));
land_mask(land_mask == 0) = nan;

params_masked = asi_regimes.param;
params_masked_land = params_masked;
for dd = 1:size(params_masked,3)
    params_now = params_masked(:,:,dd);
    params_now(land_mask==1) = -999;
    params_masked_land(:,:,dd) = params_now;
end

cm = cmocean('balance');
cm(1,:) = 0;

dd = 182;

% Add panels
ax3 = subplot(2,2,2);
pcolor(asi_regimes.lon,asi_regimes.lat,asi_regimes.param(:,:,dd)')
shading flat
cb = makecb([0 1],cm,'ASI Regime Parameter',30,ax3);
makespruce(30)
maketitle('(b) Air-Sea Interaction Regimes on July 1',40)
% Overlay the AS + BoB polygons used in panel (c)
hold on
roi_color = [0.30 0.30 0.30]; 
patch(AS_lon,  AS_lat,  'k','Edgecolor',roi_color, 'FaceColor','none', 'LineWidth',5,'linestyle','--');
patch(BoB_lon, BoB_lat, 'k','Edgecolor',roi_color, 'FaceColor','none', 'LineWidth',5,'linestyle','--');

hold off

ax4 = subplot(2,2,1);
pcolor(asi_regimes.lon,asi_regimes.lat,ws_clim_OAflux_short(:,:,dd)')
shading flat
cb = makecb([0 15],cmocean('rain'),'Windspeed (m/s)',30,ax4);
makespruce(30)
maketitle('(a) 10m Windspeed Climatology on July 1',40)

% Overall title
sgtitle('', 'interpreter', 'latex', 'fontsize', 50);
annotation('textbox', [0.2 0.95 0.6 0.04], ...
    'String', 'Northern Indian Ocean ASI Regimes vs. Windspeed', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'FontSize', 45, 'interpreter', 'latex');

% Now make correlation 

% Calculate spatial correlation
for dd = 1:size(asi_regimes.param,3)
    ws_curr = ws_clim_OAflux_short(:,:,dd);
    reg_curr = asi_regimes.param(:,:,dd);
    goodind = ~isnan(ws_curr(:)) & ~isnan(reg_curr(:));
    [cmat,pmat] = corrcoef(ws_curr(goodind),reg_curr(goodind));
    C(dd) = cmat(2,1);
    P(dd) = pmat(2,1);
end

% Plot
subplot(2,2,4)
cplot = C;
p_cutoff = BH(P,0.1); % Benjamini-Hochberg procedure in time
cplot(P>p_cutoff) = nan;
dateplot = datetime(2024,1,1):datetime(2024,1,365);
plot(dateplot,cplot,'linewidth',2)
hold on
yline(0,'linewidth',2,'linestyle','--')
ylabel('Correlation','interpreter','latex')
makespruce(30)
maketitle('(d) Windspeed-ASI Regime Spatial Correlation',40)

% Adjust x-axis ticks and format
xticks(datetime(2024,1,1):calmonths(1):datetime(2024,12,31));  % Set tick positions at the start of every month
xtickformat('MMM');  
