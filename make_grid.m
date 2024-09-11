% Function: make_grid
% This script performs spatial kriging on a grid based on cluster points.
% It generates grid occupation data and plots the results.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz  <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input Files:
% - data/seamask_norm_0.05: Seamask data
% - c14mat/C14_europe0: C14 dates and site information
% - outputDirectory/mat/clusti_<time>_120: Cluster information for each time segment

% Output Files:
% - outputDirectory/plots/grid_regnum_<dlon>_<time>.png: Plot of grid region numbers
% - outputDirectory/plots/grid<dlon>.png: Plot of the grid
% - outputDirectory/area_<dlon>_<breakPoints(1)>_<breakPoints(end)>.mat: Area data

% Variables:
% - offset: Offset for parallel computation
% - dll: Grid resolution
% - MaxOcc: Maximum occupation per grid cell
% - angl: Angles for spatial weighing
% - sn: Sign matrix for spatial weighing
% - long: Longitude grid
% - latg: Latitude grid
% - nx: Number of longitude grid points
% - ny: Number of latitude grid points
% - weigh: Spatial weighing function
% - radmax: Maximum radius for spatial weighing
% - values: Matrix for grid occupation values
% - regs: Matrix for region indices
% - clustn: Unique cluster indices
% - regionlon: Longitude of region centers
% - regionlat: Latitude of region centers
% - ncolor: Number of clusters/regions
% - ncol: Number of columns in plot
% - nrow: Number of rows in plot
% - dxp: Horizontal spacing for plot
% - dyp: Vertical spacing for plot
% - tii: Time segment index
% - ti: Current time segment
% - nreg: Number of regions


close all
clear all
load_pars; % Sets common parameters (e.g., outputDirectory, latitudeLimits, breakPoints)

% Retrieve control parameters from function argument (in parallel computation mode)
if exist('offset', 'var')
    offset = str2num(offset);
    npart = 64;
else
    offset = 1;
    npart = 1;
end

if exist('dll', 'var')
    dlon = str2num(dll);
    dlat = dlon;
else
    dlon = 0.05;
    dlat = 0.05;
end

% Set parameters and switches
MaxOcc = 4; % Maximal occupation per grid cell
angl = (0:0.04:1) * 2 * pi;
sn = [[1 1]; [1 -1]; [-1 1]; [-1 -1]];

% Grid information
long = longitudeLimits(1):dlon:longitudeLimits(end);
latg = latitudeLimits(1):dlat:latitudeLimits(end);
nx = length(long);
ny = length(latg);

% Spatial weighing function
radmax = round(1.8 / dlon);
for ix = 0:radmax-1
    for iy = 0:ix
        rad = sqrt(ix^2 + iy^2);
        if rad <= radmax + 0.01
            % Exponential decay with distance from center
            weigh(1 + ix, 1 + iy) = exp(-(7 * rad * dlon));
            weigh(1 + iy, 1 + ix) = weigh(1 + ix, 1 + iy);
        else
            weigh(1 + ix, 1 + iy) = 0;
            weigh(1 + iy, 1 + ix) = 0;
        end
    end
end
weigh(1, 1) = 1;

% Load seamask
load('data/seamask_norm_0.05'); % seamask_High

% Load C14 dates and site info
load('c14mat/C14_europe0'); % 'lonsn', 'latsn', 'C14agesn', 'C14SDsn', 'SiteIDsn', 'datIDsn'
nt = length(breakPoints);

% Number of columns and rows in plot; size
ncol = round(sqrt(nt)) + 1;
nrow = ceil(nt / ncol);
dxp = 0.97 / ncol;
dyp = 0.97 / nrow;

% Loop over time segments
tii = 0;
for ti = breakPoints
    % Clear matrices
    values = zeros(MaxOcc, nx, ny); % -9999 all ocean
    regs = zeros(MaxOcc, nx, ny); % All ocean

    % Load geo-position of sites and index that links sites to clusters
    load([outputDirectory 'mat/clusti_' num2str(ti) '_120']);
    lons = lonsn;
    lats = latsn;
    nmax = length(lons);

    % Clear index file
    clustn = unique(clusti);
    clustn = clustn(~isnan(clustn));
    clustn = clustn(clustn ~= 0);

    ncolor = length(clustn); % Number of clusters/regions
    fprintf('t=%d #regions=%d\n', ti, ncolor);

    % Center position of regions
    for i = 1:ncolor
        ii = find(clusti == i);
        ii(ii > nmax) = [];
        lon = lonsn(ii);
        lat = latsn(ii);
        regionlon(i) = mean(lon);
        regionlat(i) = mean(lat);
        fprintf('%2d:%1.1f %1.1f\t', i, regionlon(i), regionlat(i));
    end
    fprintf('\n');

    % Kriging: distribution of site positions -> grid occupation
    make_grid_regions;

    % Plot output
    figure(2);
    set(gcf, 'position', [15 20 850 750], 'Color', 'w', 'Visible', 'on');
    clf;
    subplot('Position', [0.01 0.01 0.98 0.98]);
    imagesc(flipud(reg'));
    for i = 1:ncolor
        x = (regionlon(i) - longitudeLimits(1)) / (longitudeLimits(2) - longitudeLimits(1));
        y = (regionlat(i) - latitudeLimits(1)) / (latitudeLimits(2) - latitudeLimits(1));
        annotation('textbox', 0.01 + [x y 0.05 0.05] * 0.98, 'string', [num2str(i) ':' num2str(round(area(i) * 1E-4))], ...
                   'FontWeight', 'bold', 'FontSize', 24, 'EdgeColor', 'none', 'Color', 'w');
    end
    annotation('textbox', [0.15 0.82 0.12 0.05], 'string', num2str(ti), 'FontWeight', 'bold', 'FontSize', 32, 'EdgeColor', 'none', 'Color', 'w');
    outfilename = [outputDirectory 'plots/grid_regnum_' num2str(dlon) '_' num2str(ti) '.png'];
    print('-dpng', '-r300', outfilename);

    figure(1);
    set(gcf, 'position', [10 15 1450 1050], 'Color', 'w', 'Visible', 'on');
    iy = floor(tii / ncol);
    ix = mod(tii, ncol);
    subplot('Position', [0.01 + ix * dxp 0.01 + iy * dyp dxp * 0.98 dyp * 0.98]);
    imagesc(flipud(val'));
    set(gca, 'XTick', [], 'YTick', []);
    text(100, 100, num2str(ti), 'FontWeight', 'bold', 'FontSize', 22, 'Color', 'w');
    tii = tii + 1;
    nreg(tii) = ncolor;
end

% Save areas to file
save([outputDirectory 'area_' num2str(dlon) '_' num2str(breakPoints(1)) '_' num2str(breakPoints(end)) '.mat'], 'areav', 'nreg');

% Save plot to file
set(gcf, 'PaperPositionMode', 'auto', 'InvertHardCopy', 'off', 'Visible', 'on');
outfilename = [outputDirectory 'plots/grid' num2str(dlon) '.png'];
print('-dpng', '-r300', outfilename);
