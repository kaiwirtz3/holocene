% This script loads logit model results and stores low-high-pass filtered probability
% differences into a matrix. It processes the probability of boom/bust events and
% pseudo-RGR, interpolates the data, and writes the results to a file.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input Files:
% - target_ts_0.mat: Contains 'dat' and 'legdat' variables
% - glmres_<tag>.mat: Contains logit model results for different tags

% Output Files:
% - legdat.dat: Contains the updated legend data

% Variables:
% - outputDirectory: Directory for input and output files
% - nd0: Size of old data matrix
% - vis: Cell array of model variants
% - varname: Cell array of variable names
% - nv: Index for model variants
% - tag: Model variant tag
% - sp: Sum of probabilities
% - fac: Factor for probability adjustment
% - tip1: Time vector
% - it: Indices for time vector within the range of timres
% - j: New index for extending data matrix
% - ts: Time series data for probabilities
% - tmov, toff: Parameters for moving average
% - adds: Additional description for legend
% - fid: File identifier for writing legend data

load([outputDirectory 'target_ts_0.mat']); % Load target time series data

nd0 = length(legdat); % Size of old data matrix

% Load logit model results
vis = {'-TSI_clim_area', '-TSI_clim', 'area_area', '3'};
varname = {'TSI', 'stress-free tree growth', 'climate stability', 'RGR($t$-350a)'};

for nv = 1:3
    % Name of model variant
    tag = [vis{nv}];
    load([outputDirectory 'glmres_' tag '.mat']); % Load model results

    sp = sum(prob, 2);
    fac = sp(2) / sum(sp);
    tip1 = dat(:, 1)'; % Time vector
    it = find(tip1 >= timres(1) & tip1 <= timres(end));

    % New index; extend data matrix
    j = nd0 + nv;

    % Loop over boom/bust/pseudo-RGR
    for i = 3:3 % Only pseudo-RGR
        if i <= 2
            ts = prob(i, :); % Probability of boom or bust
        else
            % Probability of boom minus probability of bust
            ts = fac * prob(1, :) - (1 - fac) * prob(2, :);
        end

        ts = interp1(timres, ts, tip1(it), 'linear', 'extrap');
        ts = movweighavg(tip1(it) * 1E3, ts, tmov, toff);

        % Refine description
        if i == 3
            dat(it, j) = ts;
            switch (nv)
                case 1
                    adds = '(climate+)';
                case 2
                    adds = '(climate)';
                case 3
                    adds = '(autoregressive lags=350a+210a)';
            end
            legdat{j} = ['logit model ' num2str(4 - nv) 'V ' adds];
        end
    end
end

% Write out variables and indices
fid = fopen('legdat.dat', 'w');
j = floor(length(legdat) / 2);

for i = 1:j
    fprintf(fid, '%2d %25s\t\t', i, legdat{i});
    fprintf(fid, '%2d %23s\n', j + i, legdat{j + i});
    fprintf('%2d %25s\t\t', i, legdat{i});
    fprintf('%2d %23s\n', j + i, legdat{j + i});
end

fclose(fid);
