% This script reads C14 dates from an Excel file, processes the data by correcting
% wrong latitude/longitude formats, removes NaNs, and saves the cleaned data into
% a MATLAB binary file. It also handles UTF-8 problems by using index vectors.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input Files:
% - c14mat/C14_dASIS.xlsx: Excel file containing C14 dates and site information

% Output Files:
% - c14mat/C14_dASIS.mat: MATLAB binary file containing processed C14 data

% Variables:
% - fname: Filename for the input Excel file
% - off: Offset for column indices in the Excel file
% - IDi: Column index for SiteID in the Excel file
% - Dbi: Column index for Dbase in the Excel file
% - lli: Column index for latitude/longitude in the Excel file
% - num, txt, raw: Data read from the Excel file
% - C14age: Vector of C14 ages
% - C14SD: Vector of C14 standard deviations
% - SiteID: Vector of site IDs
% - Dbase: Vector of database entries
% - lon, lat: Vectors of longitude and latitude
% - ii: Index vector for correcting and filtering data
% - situ, ia, ic: Unique site IDs and their indices
% - lats, lons: Processed latitude and longitude vectors
% - C14ages, C14SDs: Processed C14 age and standard deviation vectors
% - SiteIDs: Processed site ID vector
% - datIDs: Index vector for data entries

close all;

load_pars; % Sets common parameters (outputDirectory, cc, latitudeLimits, regs)

% Load own compilation (dASIS+)
fname = 'C14_dASIS';
off = 3;
IDi = 2;
Dbi = 8;
lli = 4;

% Read C14 dates from XLS file
[num, txt, raw] = xlsread(['c14mat/' fname '.xlsx']); % "labnr" "c14age" "c14std" "lat" "lon"

% Distribute data into named vectors
C14age = num(:, 1 + off);
C14SD = num(:, 2 + off);
SiteID = txt(2:end, IDi);
Dbase = txt(2:end, Dbi);
lon = num(:, lli + off);
lat = num(:, lli + 1 + off);

% Correct few wrong lon/lat value formats
ii = find(abs(lon) > 70);
if ii
    lon(ii) = lon(ii) * 1E-3;
end

ii = find(abs(lat) > 90);
if ii
    lat(ii) = lat(ii) * 1E-3;
end

% Identify and remove NaNs
ii = find(~isnan(C14age) & ~isnan(C14SD) & ~isnan(lon) & ~isnan(lat) & ~contains(Dbase, 'OUT'));
fprintf('%s: %d -> %d (no NaNs, OUT)\n', fname, length(lat), size(ii));
lat = lat(ii);
lon = lon(ii);
C14age = C14age(ii);
C14SD = C14SD(ii);
SiteID = SiteID(ii);

% Uses index vector instead of SiteID strings (because of UTF-8 problems)
SiteID = strtrim(SiteID);
[situ, ia, ic] = unique(SiteID);
SiteID = num2str(ic);
fprintf('from %d sites\n', length(ia));

lats = lat;
lons = lon;
C14ages = C14age;
C14SDs = C14SD;
SiteIDs = SiteID;
datIDs = 1:length(C14age); % Create index vector

% Save data to matlab binary
save(['c14mat/' fname], 'C14age', 'lons', 'lats', 'C14SDs', 'SiteIDs', 'datIDs');
