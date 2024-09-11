% This script processes and plots calculated growth rates (RGR) for different methods.
% It loads area information, calculates pooled and area-based RGR, smooths and detrends the data,
% and generates plots for visual comparison. The script also saves the processed data and plots to files.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input files:
% - <outputDirectory>/area_0.05_<breakPoints(1)>_<breakPoints(end)>.mat
% - <outputDirectory>/mat/AllPop_europe0<tag>_all.mat

% Output files:
% - <outputDirectory>/plots/CompRGR<tag>.png
% - <outputDirectory>/avg_rgr_all.mat

% Variables:
% - fs: Font size for plots
% - colj: Color map for plots
% - leg: Legend entries for plots
% - nregv: Number of regions
% - dt: Time step for time vectors
% - twin: Time window for linear weighing
% - ddt: Time vector for linear weighing
% - tim0, tim: Time vectors for calculations
% - rtim, tirgr: Time vectors for RGR calculations
% - wei: Weights for linear weighing
% - ntag: Number of SPD methods
% - tag, scs: Labels for SPD methods
% - rgri: Interpolated RGR for pooled data
% - rgr_arn: Interpolated RGR for area-based data
% - rgr_m: Matrix to store RGR data
% - le: Legend entries for plots
% - file: Filename for saving plots

clear all;
close all;

load_pars; % sets common parameters (e.g., outputDirectory, dataTags, breakPoints)

fs = 22;
colj = jet(6);
colj(7, :) = zeros(1, 3);
colj(8, :) = ones(1, 3) * 0.5;
colj(5, :) = [0.8 0.5 0.1];
leg = {'pooled 1yr', 'pooled w150', 'pooled smooth+detrend', 'area 1yr', 'area w150', 'area smooth+detrend'};

% load area information
load([outputDirectory 'area_0.05_' num2str(breakPoints(1)) '_' num2str(breakPoints(end)) '.mat']); % 'areav', 'nreg'
nregv = nreg;

% time vectors
dt = 25;
twin = 2 * dtb - 100;
ddt = -twin:dt:twin;
tim0 = breakPoints(1) + ddt(1);
tim = tim0:dt:(breakPoints(end) + ddt(end));
rtim = (timeLimits(1):0.01:timeLimits(2)) * 1E3;
tirgr = rtim * 1E-3;

% linear weighing in overlapping time windows
wei = -(abs(ddt) - twin) / 100;
wei(wei > 1) = 1;

ntag = 1;

% loop over SPD methods
for tagi = 1:ntag
    tag = dataTags{tagi}; % label of SPD method
    scs = tag;
    scs(regexp(scs, '[_]')) = [];

    % load pooled RGR (for Europe)
    load([outputDirectory 'mat/AllPop_europe0' tag '_all.mat']); % poptime, ym, trgr, rgr, nreg

    % bring both rgr estimates on same timeline
    rgri = interp1(flip(trgr), flip(rgr), rtim, 'linear', 'extrap');

    % calculate area based average of RGR
    calc_aravg_rgr

    % bring both rgr estimates on same timeline
    rgr_arn = interp1(tim, rgr_ar, rtim, 'linear', 'extrap');

    % open and clear figure
    gcf = figure(tagi);
    clf;
    set(gcf, 'position', [1 1 940 700], 'Color', 'w', 'Visible', 'on');

    gca = subplot('Position', [0.11 0.11 0.88 0.88]);

    % plot settings
    set(gca, 'XDir', 'reverse', 'fontsize', fs, 'Fontweight', 'bold', 'tickdir', 'out');
    set(gca, 'XLim', timeLimits, 'Box', 'on', 'YLim', [-2 5] * 1E-3, 'XTick', 3:9);
    ylabel(['Growth rate (ka' char([hex2dec('207B') hex2dec('00B9')]) ')']);
    hold on
    xlabel(['Time (kyr BP)'], 'FontName', 'Arial', 'FontSize', fs);
    text(8.5, 4.5E-3, ['RGR ' scs], 'FontSize', fs + 4, 'Fontweight', 'bold');

    % store growth rates
    nt = (tagi - 1) * 6;
    rgr_m(:, nt + 1) = rgri; % full/pooled
    rgr_m(:, nt + 4) = rgr_arn; % area based

    % smooth and detrend growth rates
    for m = 0:1
        ts = movweighavg(rtim, rgr_m(:, nt + 1 + m * 3), 152, 40); % smooth
        [ut, avgrowth1500] = movavg(tirgr, ts, 1.5);
        rgr_m(:, nt + 2 + m * 3) = ts;
        rgr_m(:, nt + 3 + m * 3) = ts - avgrowth1500;
    end

    % plot smoothed and original growth rates
    le(1) = plot(tirgr, rgr_m(:, nt + 2), '-', 'color', 'r', 'Linewidth', 3);
    le(2) = plot(tirgr, rgr_m(:, nt + 6), '-', 'color', 'k', 'Linewidth', 3);
    plot(tirgr, rgr_m(:, nt + 1), '-', 'color', 'r', 'Linewidth', 1)
    plot(tirgr, rgr_m(:, nt + 4), '-', 'color', 'k', 'Linewidth', 1)

    pl = legend(le, 'pooled', 'area based');

    % save plot to PNG file
    file = [outputDirectory 'plots/CompRGR' tag '.png'];
    set(gcf, 'PaperPositionMode', 'auto', 'Visible', 'on');
    print('-dpng', '-r300', file);
end

% save RGR time-series to file
save([outputDirectory 'avg_rgr_all'], 'tirgr', 'rgr_m', 'leg');
