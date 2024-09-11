% This script calculates the area-based average of Relative Growth Rate (RGR) over a time window.
% It loads data for each time window, processes the RGR values, and stores the weighted average
% into a vector.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Output Files:
% - None specified explicitly in the script

% Variables:
% - tii: Index for time window
% - n: Length of weight vector
% - rgr_ar: Vector to store area-based average RGR
% - breakPoints: Time window break points
% - outputDirectory: Directory for input files
% - tag: Tag for input file names
% - rgr: Relative Growth Rate values
% - trgr: Time vector for RGR
% - tarea: Total area
% - rgrm: Weighted average RGR
% - tim0: Initial time
% - ddt: Time step
% - dt: Time interval
% - i0: Initial index for time window
% - ii: Indices for time window
% - wei: Weight vector
% - areav: Cell array of area values
% - nregv: Vector of region counts
% - nregv2: Vector to store region counts for each time window

tii = 1;
n = length(wei);
rgr_ar = zeros(1, length(tim));

% Loop over time window
for ti = breakPoints

    % Load data for each time window
    load([outputDirectory 'mat/AllPop' tag '_' num2str(tii) '.mat']);
    % Variables: poptime=tm, ymv=ymv, trgr=tirgr, rgr=rgrv, nreg=nregions

    rgr(isnan(rgr)) = 0;

    if length(trgr) ~= n
        fprintf('dim mismatch %d: %d\n', length(trgr), n);
        ti = 9E9;
    else
        tarea = sum(areav{tii}(1:nregv(tii))); % Total area

        rgrm = areav{tii}(1:nregv(tii)) * rgr / tarea; % Weighted average

        % Down-values near the edges of the time window
        i0 = (ti - tim0 + ddt(1)) / dt + 1;
        ii = i0:(i0 + n - 1);
        rgr_ar(ii) = rgr_ar(ii) + fliplr(wei .* rgrm); % Stores into vector

        tii = tii + 1;
        nregv2(tii) = nreg;
    end
end
