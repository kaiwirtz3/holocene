% This script calculates the overlap of booms and busts in time series data.
% It identifies positive and negative phases/peaks, applies linear weighing,
% and computes the overlap as a weighted sum of products. The results are used
% to determine the match factor of equal periodic signals.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Variables:
% - tmaxo: Maximum time for overlap calculation
% - fi, tot, tweigh: Accumulators for overlap calculations
% - ii2: Indices for time range within tmaxo and 3
% - crit: Critical value for phase identification
% - dtt: Time step difference
% - delta: Delta value for linear weighing
% - grv: Gray value for color adjustment
% - dim: Dimension for linear weighing
% - yl: Y-limits for graphical settings
% - eps: Epsilon value for graphical settings
% - yli: Adjusted Y-limits for graphical settings
% - fa: Face alpha value for transparency
% - pcol: Color for bars
% - bgcol: Background color for bars
% - overlap: Vector to store overlap values
% - boombust: Vector to mark booms and busts
% - sgp: Sign for booms and busts
% - iy1: Index for Y-limits
% - pg, pgi: Peaks and their indices
% - jj: Number of peaks
% - xm: Time of peak
% - im: Index of peak
% - weigh: Weighing vector
% - i0, i1: Indices for start and end of phase
% - di: Dimension for linear weighing
% - i01, i10: Indices for linear weighing change
% - isp: Indices for phase
% - nl: Number of elements in phase
% - df, dft: Overlap values for weighted sum
% - tip1, ts: Time series data
% - to2: Total overlap value
% - fr: Match factor of equal periodic signals

% Coefficients for linear filter
tmaxo = 8.15;
fi = 0;
tot = 0;
tweigh = 0;
ii2 = find(tip2 <= tmaxo & tip2 >= 3);
crit = 0.5 * std(avgrde(ii2));
dtt = abs(tip2(2) - tip2(1));
delta = 0.095;
grv = 1;
dim = round(delta / dtt);

% Graphical settings
yl = get(gca, 'Ylim');
eps = abs(yl(2) - yl(1)) / 700;
yli = [yl(2) - eps eps -eps yl(1) + eps];
fa = 0.667; % 'FaceAlpha' factor

pcol = [0.6 0.2 0.]; % Bar color
bgcol = [225 200 150] / 230;

% Loop over booms and busts
for i = 1:2
    overlap = ones(1, size(dat, 1));
    boombust = zeros(1, size(dat, 1));
    sgp = (3 - 2 * i); % Sign of booms and busts
    iy1 = 2 * (i - 1) + 1;

    % Find positive or negative phases/peaks
    [pg, pgi] = findpeaks(tip2, avgrde * sgp, stdc, 9);

    % Delete peaks outside time-window
    ii2 = find(pg < 3 | pg > tmaxo);
    if ii2
        pg(ii2) = [];
        pgi(ii2) = [];
    end

    jj = length(pg);

    % Loop over phases/peaks
    for j3 = 1:jj
        xm = pg(j3); % Time
        im = pgi(j3); % Index
        weigh = zeros(1, size(dat, 1)); % Clear weighing vector

        % End of phase (smallest age kaBP)
        i0 = max(find(tip2 < xm & avgrde * sgp < crit));
        di = round(min(1, abs(avgrde(im)) / crit) * dim);
        if isempty(i0)
            i0 = im;
        end
        if i0 <= di
            i0 = 1;
        else
            i0 = i0 - di;
        end

        % Start of phase
        i1 = min(find(tip2 > xm & avgrde * sgp < crit));
        if isempty(i1)
            i1 = im;
        end
        if i1 >= length(tip2) - di
            i1 = length(tip2);
        else
            i1 = i1 + di;
        end

        % Linear weighing change at edge
        weigh(i0:(im - 1)) = ((i0:(im - 1)) - i0) / dim;
        weigh(im:i1) = (i1 - (im:i1)) / dim;
        i01 = min(find(weigh >= 1));
        if isempty(i01)
            i01 = im;
        end
        i10 = max(find(weigh >= 1));
        if isempty(i10)
            i10 = im;
        end
        weigh(weigh > 1) = 1;

        tweigh = tweigh + sum(weigh);
        isp = i0:i1;
        isp = isp'; % Vector of indices
        nl = length(isp);

        boombust(min((isp)):max((isp))) = 1; % Mark phases

        % Add background bars
        add_bars

        % Overlap: weighed sum of product
        df = nansum(abs(tsi .* weigh(isp1) .* overlap(isp1)));
        fi = fi + df;
        dft = nansum(abs(tsi .* overlap(isp1)));
        tot = tot + dft;

        overlap(isp) = 0; % Avoid double counting

    end % for j3

end % for i=1:2

ii2 = find(tip1 < tmaxo & tip1 >= 3);
to2 = nansum(abs(ts(ii2)));
fr = fi / (to2 * 0.92); % Match factor of equal periodic signals
