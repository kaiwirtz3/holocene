% This script performs kriging to distribute site positions into grid occupations.
% It processes clusters/regions, checks for singularity points, calculates grid index positions,
% and normalizes value fields. The script also calculates the area of each region and stores
% region information in a binary file.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz  <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Output Files:
% - regiongrid_<ti>.mat: Binary file containing region information

% Variables:
% - vc: Minimum value for cluster contribution
% - minc: Minimum distance between sites in the same cluster
% - ncolor: Number of clusters/regions
% - ii: Indices of sites in the current cluster
% - ip: Index of the current site in the cluster
% - mind: Minimum distance to the nearest neighbor site
% - jp: Index for nearest neighbor search
% - ix, iy: Grid index positions
% - iso: Indicator for grid cell occupation
% - nv: Occupation layer index
% - dd: Directional index
% - sea: Indicator for marine straits
% - occ: Occupation indicator
% - jxi, jyi: Indices for neighboring grid cells
% - jx, jy: Offsets for neighboring grid cells
% - x1, y1: Coordinates of neighboring grid cells
% - mv, mi: Minimum value and its index in the occupation layers
% - ind1, ind: Indices for valid and invalid grid cells
% - vmax: Maximum value for normalization
% - sval: Sum of values for normalization
% - val, reg: Squeezed value and region matrices
% - arf: Area factor for region calculation
% - area: Array of region areas
% - areav: Cell array to store area information

vc = 1;
minc = 100;

% Loop over clusters/regions
for i = 1:ncolor
    if exist('clustdat')
        ii = find(clusti == clustn(clustdat(i, 2)));
    else
        ii = find(clusti == i);
    end
    ii(find(ii > length(lats))) = [];

    % Loop over sites in cluster
    for ip = 1:length(ii)

        % Check for singularity points
        % After clustering, few sites can be located within or close to 'neighbor' clusters

        % Distance to nearest neighbor site of the same cluster should be minc (=100km)
        mind = 9E9;
        jp = 1;

        % Only screen north and east regions for singularity points
        if (regionlon(i) > 30 || regionlat(i) > 59)
            while mind > minc && jp <= length(ii)
                cd = cl_distance(lons(ii(ip)), lats(ii(ip)), lons(ii(jp)), lats(ii(jp)));
                if cd < mind && cd > 2 && ip ~= jp
                    mind = cd;
                end
                jp = jp + 1;
            end
        end

        if jp <= length(ii)

            % Calculate grid index position from site geo-location
            ix = 1 + floor((lons(ii(ip)) - long(1)) / dlon);
            iy = 1 + floor((lats(ii(ip)) - latg(1)) / dlat);

            % Valid grid index?
            if ix <= nx && iy <= ny && ix > 0 && iy > 0
                iso = 0;

                % Occupation of grid cell by other clusters
                for nv = 1:MaxOcc
                    iso = iso + (regs(nv, ix, iy) ~= i || values(nv, ix, iy) < vc);
                end

                % Full occupation by others or small contribution by own cluster
                if iso == 4 && value(ix, iy) > -1E-3
                    for dd = 1:4
                        sea = 0;
                        occ = 0;

                        % Check for marine straits
                        for jxi = 1 + (sn(dd, 1) < 0):radmax * 0.7
                            jx = jxi - 1;
                            for jyi = 1 + (sn(dd, 2) < 0):radmax * 0.7
                                jy = jyi - 1;
                                x1 = ix + jx * sn(dd, 1);
                                y1 = iy + jy * sn(dd, 2);
                                if x1 <= nx && x1 > 0 && y1 <= ny && y1 > 0
                                    if value(x1, y1) < 0
                                        sea = sea + 1 / (jx * jx + jy * jy);
                                    end
                                end
                            end
                        end

                        % Loop over neighboring grid cells - x
                        for jxi = 1 + (sn(dd, 1) < 0):radmax
                            jx = jxi - 1;

                            % Loop over neighboring grid cells - y
                            for jyi = 1 + (sn(dd, 2) < 0):radmax
                                jy = jyi - 1;
                                if weigh(jxi, jyi) / (sea * sea * sea + 1) > 5E-4 % Inside radius
                                    x1 = ix + jx * sn(dd, 1);
                                    y1 = iy + jy * sn(dd, 2);
                                    if x1 <= nx && x1 > 0 && y1 <= ny && y1 > 0 % Check for domain
                                        if (value(x1, y1) >= 0)
                                            iso = 1;
                                            nv = 1;
                                            while iso == 1 && nv <= MaxOcc
                                                % Occupied cell by 1st
                                                iso = (regs(nv, x1, y1) ~= i && regs(nv, x1, y1) > 0);
                                                nv = nv + 1;
                                            end
                                            if iso == 0 && nv == MaxOcc + 1 % All 4 layers occupied: find lowest layer
                                                [mv, mi] = min(values(:, x1, y1));
                                                if mv < weigh(jxi, jyi)
                                                    regs(mi, x1, y1) = i;
                                                    values(mi, x1, y1) = weigh(jxi, jyi);
                                                end
                                            else % Set a new entry
                                                nv = nv - 1;
                                                values(nv, x1, y1) = values(nv, x1, y1) * (regs(nv, x1, y1) == i) + weigh(jxi, jyi);
                                                regs(nv, x1, y1) = i;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                else
                    fprintf('rm %d %d\t %1.1f/%1.1f\t%1.1f/%1.1f\t%1.3f\n', i, ip, lons(ii(ip)), regionlon(i), lats(ii(ip)), regionlat(i), mind);
                end
            end
        end
    end
end

for i = 0:-ncolor
    fprintf('%d %d %d\n', i, length(find(regs(1, :, :) == i)), length(find(regs(2, :, :) == i)));
end

ind1 = find(value >= 0);
vmax = 4;
for nv = 1:MaxOcc
    ind = find(values(nv, :, :) > vmax);
    values(nv, ind) = vmax;
end

% Normalize all value fields
sval = (sum(values, 1)) + 1E-5;
ind = find(value < 0);
for nv = 1:MaxOcc
    values(nv, :, :) = values(nv, :, :) ./ sval;
    % Mask sea areas
    values(nv, ind) = NaN;
    regs(nv, ind) = 0;
end
fprintf('MEAN %1.3f\tMAX %1.3f\n', nanmean(nanmean(values(1, ind1))), max(max(values(1, ind1))));

val = squeeze(values(1, :, :));
reg = squeeze(regs(1, :, :));
ind = find(val < 0);
val(ind) = NaN;
reg(ind) = 0;

% Sort: bring higher weight to front
for nv = 2:MaxOcc
    ind = find(squeeze(values(nv, :, :)) > val);
    if ind
        reg(ind) = squeeze(regs(nv, ind));
        val(ind) = values(nv, ind);
    end
end

% Calculate area of each region
clear area
for i = 1:ncolor
    x = regionlon(i);
    y = regionlat(i);
    arf = cl_distance(x, y, x, y + dlat);
    arf = arf * cl_distance(x, y, x + dlon, y);
    ind = find(reg == i);
    area(i) = arf * length(ind);
    fprintf('%2d %1.0f\t%1.2f * %d\n', i, area(i), arf, length(ind));
end
areav{tii + 1} = area;
