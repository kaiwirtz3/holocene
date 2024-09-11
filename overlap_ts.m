% This script calculates and plots the overlap between time-series data.
% It processes data from logit model results, computes statistical measures,
% and generates plots for different RGR reconstructions and compared time-series.
% The script also saves the processed data and plots to files.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input files:
% - legdat.dat (or equivalent data source)

% Output files:
% - <outputDirectory>/plots/RGR_<str>_<stl>.png
% - <outputDirectory>/target_ts_<length(spv)>.mat

% Variables:
% - timeLimits: Time range for analysis
% - tmax: Maximum time for statistical measures and graphical output
% - tavg: Time window for moving average
% - fs: Font size for plots
% - y0: Y-offset for annotations
% - sdf: Scaling factor for the second curve
% - yl: Y-axis limits
% - col2: Color for the second plot
% - spv, spr: Indices of compared time-series and RGR reconstructions
% - stdc: Standard deviation constant
% - nrgr: Number of RGR reconstructions
% - dt: Time step for common time vector
% - time: Common time vector
% - tip1, tip2: Time points for data
% - avgrde: Average growth rate data
% - legdat: Legend data
% - gcf, gca: Figure and axis handles for plotting
% - xl, xlp: X-label and its position
% - sc: Scaling factors
% - j, cvar: Indices and variables for comparison
% - ts: Time-series data
% - le: Legend entries
% - leg_str: Legend strings
% - stl: String for legend
% - r, p: Correlation coefficient and p-value
% - stat1, statnam: Statistical measures and names
% - file: Filename for saving data and plots

load_pars; % sets common parameters (outputDirectory, cc, latitudeLimits, regs)
timeLimits = [2.8 10.2];
tmax = 8.9; % max time for statistical measures (ka BP) and graphical output
tavg = 1.5; % time window moving average

% load logit model results and store low-high-pass filtered probability difference to data matrix
add_logitres

% graphical parameters
fs = 22;          % fontsize
y0 = 0.18;        % y-offset
sdf = 0.7;        % scaling of 2nd curve
yl = [-1.4 1.85]; % limits
col2 = [0.95 0.4 0.1]; % colour of 2nd plot

spv = [6:18 ones(1,6)*10]; % indices of compared TS, see legdat array or txt file 'legdat.dat'
spr = [ones(1,13)*2 4:9]; % indices of RGR reconstr., see legdat
stdc = 9E-5;
nrgr = length(spr);

% common time vector
dt = 0.010;
time = timeLimits(1):dt:timeLimits(2);
tip1 = dat(:,1)';

% loop over RGR reconstructions and compared TS
for ic = 1:nrgr
    gcf = figure(ic);
    set(gcf, 'position', [1 ic*25 1020 500], 'Color', 'w', 'Visible', 'on');
    clf;

    % mean growth rate
    ioff = spr(ic);
    avgrde = dat(:, ioff);
    if nanstd(avgrde) > 0.1
        avgrde = avgrde * 1E-3; % change units to ka-1
    end
    ii = find(~isnan(avgrde') & tip1 <= timeLimits(2) & tip1 >= 3);
    tip2 = tip1(ii)';
    avgrde = avgrde(ii);
    str = regexprep(legdat{ioff}, '_', '');

    % plot settings
    gca = subplot('Position', [0.1 0.14 0.82 0.85]);
    set(gca, 'XDir', 'reverse', 'fontsize', fs, 'Fontweight', 'bold', 'tickdir', 'out');
    set(gca, 'XLim', [3 tmax], 'Box', 'on', 'YLim', yl, 'YTick', -4:5, 'XTick', 3:9);
    if ioff < 6
        ylabel(['Growth rate (ka' char([hex2dec('207B') hex2dec('00B9')]) ')']);
    end
    hold on
    xl = xlabel(['Time (ka BP)'], 'FontName', 'Arial', 'FontSize', fs);
    xlp = get(xl, 'Position');
    xlp(2) = xlp(2) + 0.042;
    set(xl, 'Position', xlp);

    % add zero-line
    plot(timeLimits, zeros(2, 1), '-', 'Color', ones(3, 1) * 0.5, 'LineWidth', 1);

    % add text annotation
    text(6.3, 2.1, str, 'fontsize', fs, 'Fontweight', 'bold');

    sc = zeros(2, 1);
    j = spv(ic);
    cvar = dat(ii, j);
    scal = sdf / nanstd(cvar);
    ts = cvar * scal;

    % calculate overlap
    calc_overlap % calculate overlap and mark phases by bars

    % plot RGR
    le(1) = plot(tip2, avgrde * 1E3, '-', 'color', 'k', 'Linewidth', 3);
    if ioff < 6
        leg_str{1} = ['RGR ' str];
        if strfind(str, 'area')
            leg_str{1} = 'RGR Europe';
        end
    else
        leg_str{1} = str;
    end
    if ioff == -5 % add RGR Europe for other RGRs
        le(2) = plot(tip2, dat(ii, 2) * 1E3, '-', 'Color', ones(3, 1) * 0.66, 'LineWidth', 2);
        j20 = 1;
        leg_str{2} = 'RGR Europe';
    else
        j20 = 0;
    end

    % prepare 2nd yaxis
    yyaxis right
    set(gca, 'YColor', col2, 'YLim', yl);
    set(gca, 'YTick', 2 * [-1 0 1] * sdf, 'YTicklabel', ['-2'; '0 '; '+2'], 'fontsize', fs, 'FontWeight', 'b');
    yla = ylabel('normalized', 'fontsize', fs, 'FontWeight', 'b');
    yla.Position = yla.Position + [0.05 -0.2 0];

    % plot 2nd variable such as climate variability
    stl = regexprep(legdat{j}, '_', '');
    le(2 + j20) = plot(tip2, cvar * scal, '-', 'color', col2, 'Linewidth', 5);
    le(2 + j20).Color(4) = 0.75;

    % calculate (trimmed) correlation
    x1 = avgrde;
    x2 = cvar;
    ii1 = find(~isnan(x1) & ~isnan(x2) & dat(ii, 1) <= tmax & dat(ii, 1) >= 3);
    [r, p1] = corrcoef(x1(ii1), x2(ii1));
    r = r(1, 2);
    p = p1(1, 2);

    % store stats
    stat1(ic, :) = [r * r fr * 100];
    statnam{ic} = [str ' : ' stl];
    leg_str{2 + j20} = stl;
    annotation('textbox', [0.823 y0 + 0.011 0.05 0.022], 'String', [num2str(fr * 100, '%2.0f') '%'], 'FontName', 'Arial', 'FontSize', fs - 2, 'Fontweight', 'b', 'color', col2, 'EdgeColor', 'none');

    % add and re-position legend
    pl = legend(le, leg_str);
    ny = length(le);
    set(pl, 'box', 'off', 'fontSize', fs, 'position', [0.28 0.9 0.62 0.08], 'Orientation', 'Horizontal');
    stl = regexprep(stl, ' ', '');

    % save plot to PNG
    file = sprintf('%splots/RGR_%s_%s.png', outputDirectory, str, stl);
    file = regexprep(file, ' ', '');
    set(gcf, 'PaperPositionMode', 'auto', 'InvertHardCopy', 'off');
    print('-dpng', '-r300', file);

end % for ic loop over RGR recon

% save data and related info to MATLAB binary
file = sprintf('%starget_ts_%d.mat', outputDirectory, length(spv));
fprintf('save all relevant data in %s\n', file)
save('-v6', file, 'dat', 'legdat', 'stat1', 'statnam', 'tmax');
