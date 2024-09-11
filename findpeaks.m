% This function identifies peaks in a given time series data based on a specified threshold.
% It checks for correct values and signs, finds distant peaks, and returns the peak times
% and their corresponding indices.
%
% SPDX-FileCopyrightText: 2021-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Carsten Lemmen <carsten.lemmen@hereon.de>
% SPDX-License-Identifier: CC0-1.0

% Variables:
% - tim: Time vector
% - val: Value vector
% - pc: Threshold for peak detection
% - tm: Time margin for peak detection
% - ii: Indices of values exceeding the threshold
% - di: Differences between consecutive indices
% - inext: Indices of distant peaks
% - i0: Initial index for peak detection
% - ptimes: Vector to store peak times
% - index: Vector to store indices of peaks
% - t0, t1: Start and end times for peak detection
% - ini: Indices for current peak detection
% - nini: Number of elements in ini
% - m, in: Maximum value and its index in ini
% - ini2: Indices for secondary peak detection

function [ptimes, index] = findpeaks(tim, val, pc, tm)

if ~exist('tm', 'var')
    tm = 7.6;
end

% Check for correct value and sign
ii = find((abs(val) > (pc)) & sign(val) == 1);
i = 1;
ptimes = [];
index = [];

if ii
    di = ii(2:end) - ii(1:end-1);
    inext = find(di > 5); % Distant
    i0 = 1;

    for i = 1:length(inext)
        t0 = tim(ii(i0));
        t1 = tim(ii(inext(i)));
        ini = ii(i0):ii(inext(i));
        nini = length(ini);

        [m, in] = max(abs(val(ini)));
        index = [index ini(in)];
        ptimes = [ptimes tim(ini(in))];

        ini2 = [];

        if in <= 0.39 * nini
            ini2 = ini(round(nini * 0.55):end);
        end

        if in >= 0.61 * nini
            ini2 = ini(1:round(nini * 0.45));
        end

        if ini2
            [m, in] = max(abs(val(ini2)));
            index = [index ini2(in)];
            ptimes = [ptimes tim(ini2(in))];
        end

        i0 = inext(i) + 1;
    end

    t0 = tim(ii(i0));
    t1 = tim(ii(end));
    ptimes = [ptimes mean([t0 t1])];
    index = [index round(mean([ii(i0), ii(end)]))];
end

return;
end
