% This script merges Dynamic Time Warping (DTW) time segments data and includes other time series
% such as Relative Growth Rate (RGR) and solar forcing. It processes and re-grids various datasets
% onto a common time vector, applies smoothing and detrending, and merges DTW and PCA results.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input Files:
% - avg_rgr_all.mat: Contains RGR data for Europe
% - AllPop_<continent>_NoNorm_Bin100_all.mat: Contains RGR data for continents
% - SA_spd_rgr.mat: Contains RGR data for South America
% - Steinhilber2012_Solar.dat: Contains Total Solar Irradiance data
% - bog_std.mat: Contains Northern Irish bog data
% - dtwpca2_<tmax>_<dtw_Dist_crit>_<tol*100>_<ln>.mat: Contains DTW and PCA results

% Output Files:
% - dtwpca_proxydata_<tmax>_<dtw_Dist_crit>_<tol*100>_<ln>.mat: Contains merged DTW and PCA results
% - pca_0.mat: Contains PCA data
% - target_ts_0.mat: Contains combined data

% Variables:
% - timeLimits: Time limits for the analysis
% - tmov, toff: Smoothing parameters
% - tavg: Detrend parameter
% - dt, di: Time step and interval
% - time: Common time vector
% - dat: Data matrix
% - legdat: Legend data
% - ri: Indices for RGR methods
% - contname: Names of continents
% - cc: Current continent name
% - file: File name for input data
% - trgr, rgr: Time and RGR data
% - it: Indices for time vector
% - ts: Time series data
% - ts_time, ts_m: Time and TSI data
% - ut, avg1500: Variables for moving average
% - ti, ringw_sm, bogstd_space, bogstd_time, flood_m: Bog data variables
% - split, dtw_Dist_crit, tol: DTW and PCA parameters
% - datt, datm, datall: Data matrices for DTW and PCA results
% - datm1, datall1, datm2, datall2: Temporary data matrices
% - pca, datma: PCA data matrices
% - tmov2: Smoothing parameter for PCA
% - ln: Index for RGR method
% - tj: Time segment index
% - tmax: Maximum time for DTW and PCA
% - exp_vara: Explained variance for PCA
% - t_change, pca_change: Time and PCA change data
% - pca_cm: PCA change data after smoothing
% - time20, pcats: Time and PCA data
% - itv: Indices for time segments
% - ii: Indices for intersecting time segments
% - nt: Number of time points in intersecting segments
% - ff: Factor for merging time segments
% - clim_stability: Climate stability data

clear all;
close all;

load_pars; % Sets common parameters (outputDirectory, cc, latitudeLimits, regs)
timeLimits = [2.8 10.2];

% Parameters for band-pass filtering
tmov = 151;
toff = 40; % Smoothing parameters
tavg = 1.5; % Detrend

% Common time vector and data matrix
dt = 0.01;
di = 3;
time = timeLimits(1):dt:timeLimits(2);
dat = zeros(length(time), 12) + NaN;
dat(:, 1) = time;
legdat{1} = 'time';

% RGR Europe (methods)
load([outputDirectory 'avg_rgr_all']); % Variables: leg rgr_m tirgr
ri = [6]; % 6 area-based

for i = 1:length(ri)
    j = ri(i);
    legdat{1 + i} = leg{j};

    % Re-grid on common time vector
    it = 1:min(find(time >= tirgr(end)));
    ts = interp1(tirgr, rgr_m(:, j), time(it), 'linear', 'extrap');
    ts = movweighavg(time(it) * 1E3, ts, tmov, toff);
    dat(it, 1 + i) = ts;
end

legdat{3} = 'void';
i0 = 4;

% RGR of continents
contname = {'EAsia', 'NAmerica', 'SAmerica', 'Africa', 'Australia'};

for i = 1:length(contname)
    cc = (contname{i});
    file = [outputDirectory 'mat/AllPop_' cc '_NoNorm_Bin100_all.mat'];

    if exist(file)
        load(file); % Variables: poptime=tm, ymv=ymv, trgr=tirgr, rgr=rgrv, nreg=nregions
        trgr = trgr * 1E-3;
        trgr = flipud(trgr);
        rgr = flipud(rgr);

        % Bring both rgr estimates on same timeline
        it = find(time >= (trgr(1)) & time <= (trgr(end)));
        ts = interp1(trgr, rgr, time(it), 'linear', 'extrap');
        ts = movweighavg(time(it) * 1E3, ts, tmov, toff) * 1E3;
        dat(it, i0) = movweighavg(time(it) * 1E3, ts, tmov, toff); % Second smooth
    end

    legdat{i0} = ['RGR ' cc];
    i0 = i0 + 1;
end

% RGR South America
load(['data/SA_spd_rgr']); % Variables: sa_rtim, spd1, sa_rgr
it = find(time >= sa_rtim(1) & time <= sa_rtim(end));
ts = interp1(sa_rtim, sa_rgr, time(it), 'linear', 'extrap');
ts = movweighavg(time(it) * 1E3, ts, tmov, toff);
dat(it, i0) = ts;
legdat{i0} = 'RGR South America';
i0 = i0 + 1;

% Total Solar Irradiance (TSI) from Steinhilber et al. 2012
ts = load(['data/Steinhilber2012_Solar.dat']);
ts_time = ts(:, 1);
it = 1:min(find(time >= ts_time(end)));
ts_m = interp1(ts_time, ts(:, 2), time(it), 'linear', 'extrap');
[ut, avg1500] = movavg(time(it), ts_m, 1.);
ts_m = ts_m - avg1500;
ts_m = ts_m / nanstd(ts_m(find(time(it) < 8.2)));
dat(it, i0) = -movweighavg(time(it) * 1E3, ts_m, tmov, toff); % Revert in sign!
legdat{i0} = '- TSI';
i0 = i0 + 1;

% Northern Irish bog data provided by Rowan McLaughlin
load(['data/bog_std']); % Variables: ti, ringw_sm; bogstd_space; bogstd_time;
it = find(time <= ti(end) & time >= ti(1));
ts = interp1(ti, ringw_sm, time(it), 'linear', 'extrap');
ts0 = interp1(ti, bogstd_space, time(it), 'linear', 'extrap');
ts1 = interp1(ti, bogstd_time, time(it), 'linear', 'extrap');
ts2 = interp1(ti, flood_m, time(it), 'linear', 'extrap');
dat(it, i0) = (ts - mean(ts)) / std(ts);
dat(it, i0 + 1) = -(ts0 - mean(ts0)) / std(ts0);
legdat{i0} = ['bog tree ring width'];
legdat{i0 + 1} = ['tree growth homogeneity'];
i0 = i0 + 2;

% Merge DTW & PCA results for two time-segments
ri = [7];
jj = 1;
split = 2;
dtw_Dist_crit = 90;
tol = 0.1;

% Clear all fields
datt = zeros(split, length(time), length(ri)) + NaN;
datm = zeros(split, length(time)) + NaN;
datall = zeros(1, length(time)) + NaN;
datm1 = datm;
datall1 = datall;
datm2 = datm;
datall2 = datall;
pca = zeros(5, length(time)) + NaN;
datma = zeros(5, split, length(time)) + NaN;
tmov2 = 181;
ln = ri(jj);

% Loop over time-segments
tj = 1;

for tmax = [6.2 9.5]
    file = sprintf('%sdtwpca/dtwpca2_%3.2f_%2.0f_%1.0f_%d.mat', outputDirectory, tmax, dtw_Dist_crit, tol * 100, ln);
    fprintf('tm=%1.1f %d\t loading dtw-pca from %s\n', tmax, ln, file);

    if exist(file)
        load(file);
        exp_vara(tj,:) = exp_var;
        it = max(find(time <= (t_change(1)))):min(find(time >= (t_change(end))));

        pca_cm = movweighavg(t_change * 1E3, pca_change, tmov, toff);
        datt(tj, it, j) = interp1(t_change, pca_cm, time(it), 'linear', 'extrap');

        for jp = 1:5
            pcat = movweighavg(time20 * 1E3, pcats{jp}, tmov, toff);
            datma(jp, tj, it) = interp1(time20, pcat, time(it), 'linear', 'extrap');
        end

        datm(tj, it) = nanmean(datt(tj, it, :), 3);
        itv{tj} = it;
    end

    tj = tj + 1;
end

% Second loop over time-segments
for tj = 1:2
    datall(itv{tj}) = datm(tj, itv{tj});

    for jp = 1:5
        pca(jp, itv{tj}) = datma(jp, tj, itv{tj});
    end
end

ii = intersect(itv{1}, itv{2});
nt = length(ii);
ff = (ii(end) - ii) / (nt - 1);
datm(1, ii(end)) = 0;
datm(2, ii(1)) = 0;
datall(ii) = ff .* datm(1, ii) + (1 - ff) .* datm(2, ii);

for jp = 1:5
    datma(jp, 1, ii(end)) = 0;
    datma(jp, 2, ii(1)) = 0;
    pca(jp, ii) = ff .* squeeze(datma(jp, 1, ii))' + (1 - ff) .* squeeze(datma(jp, 2, ii))';
end

[ut, avg1500] = movavg(time, datall, 1);
datall = datall - avg1500;
clim_stability = datall / nanstd(datall);
dat(:, i0) = clim_stability;

for jp = 1:2
    dat(:, i0 + jp) = pca(jp, :) / nanstd(pca(jp, :));
    legdat{i0 + jp} = ['climate PCA' num2str(jp)];
end

legdat{i0} = 'climate stability';
file = sprintf('%sdtwpca/dtwpca_proxydata_%3.2f_%2.0f_%1.0f_%d.mat', outputDirectory, tmax, dtw_Dist_crit, tol * 100, ln);
file = sprintf('%spca_0.mat', outputDirectory);
fprintf('save PCs/dPCs data in %s\n', file);
save('-v6',file,'time','pca','exp_vara','clim_stability');

% --------------------------------------
% save combined data to file
file = sprintf('%starget_ts_0.mat', outputDirectory);
fprintf('save all data in %s\n',file)
save('-v6',file,'dat','legdat','tmov','toff');
