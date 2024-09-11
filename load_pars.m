% This script sets common parameters for the study.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz  <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Add path to the m_map toolbox
addpath('~/tools/m_map');

% Geographic boundary of the study area
latitudeLimits = [34 71]; % Entire Europe
longitudeLimits = [-12 37];

% Time period of interest in kyrBP
timeLimits = [3 9.3];
timeStep = 200;
breakPoints = 3000:(2*timeStep):9800;

% Names of data treatments in SPD generation
dataTags = {'_NoNorm_Bin100', '_Norm_Bin100'};

% Output directory
outputDirectory = 'out/';

% Work directory
currentDirectory = pwd;
workDirectory = [currentDirectory '/'];

% Color settings
color1 = [223 233 130]/256;
markerColors(1,:) = [1 0.5 0];
markerColors(2,:) = [0 0.5 1];
alphabet = 'alphabetd';
