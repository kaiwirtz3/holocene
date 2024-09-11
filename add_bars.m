% This script adds background bars to a plot, representing different phases and overlaps
% with a second variable. It processes the data to create patches for left, right, and central
% bars, and handles the overlap of phases with another variable.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Variables:
% - xx: Reshaped time points for phases
% - ii2: Indices for edges
% - nli: Number of elements for edges
% - x, y: Coordinates for patches
% - ii20: Last index for left edge
% - i3: Loop index for edges
% - col: Color for patches
% - fac: Factor for transparency adjustment
% - isp1: Indices for overlap with second variable
% - tsi: Time series data for overlap
% - ix: Indices for sign check
% - nl1: Number of elements for overlap
% - iic: Indices for central bar

% ------- RGR phases
xx = tip2(isp);
xx = reshape(xx, 1, length(xx));

% ------- Left edge
ii2 = (i0:i01) - (i0 - 2);
nli = length(ii2);
x = [xx(ii2) fliplr(xx(ii2))];
y = [yli(iy1) * ones(1, 2) yli(iy1 + 1) * ones(1, 2)];
ii20 = ii2(end);

for i3 = 1:nli-1
    col = bgcol * grv;
    fac = (i3 - 1) / (nli - 1);
    patch('XData', [xx(ii2(i3:(i3 + 1))) fliplr(xx(ii2(i3:(i3 + 1))))], ...
          'YData', y, 'Facecolor', col, 'EdgeColor', 'none', 'FaceAlpha', fa * fac);
end

% ------- Right edge
ii2 = (i10:i1) - (i1 - nl + 1);
ii2(ii2 < ii20) = [];
ii2(ii2 <= 0) = [];
nli = length(ii2);

for i3 = 1:nli-1
    col = bgcol * grv;
    fac = (nli - 1 - i3) / (nli - 1);
    patch('XData', [xx(ii2(i3:(i3 + 1))) fliplr(xx(ii2(i3:(i3 + 1))))], ...
          'YData', y, 'Facecolor', col, 'EdgeColor', 'none', 'FaceAlpha', fa * fac);
end

% ------- Central bar
if (i10 > i01 + 1)
    ii2 = ((i01 + 1):(i10 - 1)) - (i0 - 1);
    ii2(ii2 <= 0) = [];
    nli = length(ii2);
    patch([xx(ii2) fliplr(xx(ii2))], ...
          [yli(iy1) * ones(1, nli) yli(iy1 + 1) * ones(1, nli)], '-', ...
          'Facecolor', bgcol * grv, 'EdgeColor', 'none', 'FaceAlpha', fa);
end

% ------- Overlap of RGR phases with 2nd variable
isp1 = find(tip2 >= min(tip2(isp)) & tip2 <= max(tip2(isp)));
ts = reshape(ts, 1, length(ts));
tsi = ts(isp1);

ix = find(sign(tsi) ~= sgp);
tsi(ix) = 0;
xx = tip2(isp1);
xx = reshape(xx, 1, length(xx));
nl1 = length(isp1);

iic = 4:(length(xx) - 3);
nl1 = length(iic);
patch([xx(iic) fliplr(xx(iic))], [tsi(iic) zeros(1, nl1)], '-', ...
      'Facecolor', pcol, 'EdgeColor', 'none', 'FaceAlpha', 0.1);
