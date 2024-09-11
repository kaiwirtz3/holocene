% This script calculates and plots the Fast Fourier Transform (FFT) of time-series data.
% It processes data from a MATLAB file, computes the FFT, and generates plots for different variables.
% The script also saves the processed data and plots to files.
%
% SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
% SPDX-License-Identifier: GPL-3.0-or-later

% Input files:
% - <outputDirectory>/target_ts_19.mat

% Output files:
% - <outputDirectory>/plots/fft_<subplotIndex>.png

% Variables:
% - labelLetters: Letters for subplot labels
% - numColumns, numRows: Number of columns and rows for subplots
% - plotWidth, plotHeight: Width and height of each subplot
% - lineColors, lineColors2: Colors for plotting lines
% - fontSize: Font size for plots
% - tickMarks: Tick marks for the x-axis
% - timeMax: Maximum time for valid data
% - selectedVariables: Indices of selected variables for FFT
% - timePeriod: Time periods for FFT calculation
% - figureHandle: Handle for the figure
% - subplotIndex: Index for the current subplot
% - colIndex, rowIndex: Column and row indices for subplots
% - legendHandles, legendLabels: Handles and labels for the legend
% - subplotPosition: Position of the subplot
% - validIndices: Indices of valid data points
% - numSubVars: Number of sub-variables
% - currentVar: Current variable index
% - dataIndices: Indices of data points for the current variable
% - timeSeriesData: Time-series data for the current variable
% - labelStr: Label string for the current variable
% - fixedPlot: Flag for fixed plot
% - plotColor: Color for the plot
% - samplingFreq: Sampling frequency for FFT
% - numPoints: Number of data points
% - fftData: FFT data
% - powerDensity: Power density of the FFT
% - freq: Frequency vector
% - logPowerDensity: Logarithm of the power density
% - period: Period vector
% - power: Power vector
% - interpolatedPower: Interpolated power for plotting
% - smoothedPower: Smoothed power for plotting
% - timePeriod0, powerSpectrum0: Time period and power spectrum for reference
% - legendHandle: Handle for the legend
% - outputFile: Filename for saving the plot

% Load common parameters (outputDirectory, ...)
load_pars;

% Load time-series data from MATLAB file
dataFile = sprintf('%starget_ts_19.mat', outputDirectory);
load(dataFile); % 'dat', 'legdat'

% Graphical settings
labelLetters = 'alphabetdef';
numColumns = 2;
numRows = 3;
plotWidth = 0.86 / numColumns;
plotHeight = 0.87 / numRows;
lineColors = [
    [0.6 0.3 0.35]; [0 0.2 0.9]; [0.3 0.3 0.9];
    [0.8 0 0.4]; [0 0.8 0.4]
];
lineColors2 = [
    [0.9 0.6 0.25]; [0.65 0 0.3]; [0 0 0];
    [0.7 0.1 1]; [0.2 0.7 0.3]; [0.1 0.4 0.8];
    [0.2 0.52 0.95]
];

fontSize = 24;
tickMarks = [260 360 500 680 950];
timeMax = 8.8;

selectedVariables = [
    [19 18 2]; [17 20 22]; [24 23 34];
    [33 29 -1]; [16 35 -1]; [30 31 -1]
];  % ATTENTION: make sure to hit the right indices of the input data matrix

% Increase resolution at larger periodicity
timePeriod = [140:10:400 430:20:830 900:40:1500];

% Create figure
figureHandle = figure(1);
clf;
set(figureHandle, 'position', [2 1 920 1020], 'Color', 'w', 'Visible', 'on');

for subplotIndex = 0:(size(selectedVariables, 1) - 1)
    colIndex = mod(subplotIndex, numColumns);
    rowIndex = numRows - 1 - floor(subplotIndex / numColumns);
    clear legendHandles;
    clear legendLabels;

    % Create subplot
    subplotPosition = [
        0.08 + colIndex * 1.1 * plotWidth,
        0.08 + rowIndex * 1.05 * plotHeight,
        plotWidth,
        plotHeight
    ];

    gca = subplot('Position', subplotPosition);
    set(gca, 'XTick', tickMarks, 'XLim', [195 1250], 'YLim', [-2.65 0.1], ...
        'tickdir', 'out', 'XScale', 'log', 'YTick', [-3:1:0], ...
        'fontsize', fontSize, 'Fontweight', 'b', 'Box', 'on');
    hold on;

    if rowIndex == 0
        xlabel('Period (a)', 'FontName', 'Arial', 'FontSize', fontSize);
    else
        set(gca, 'XTickLabel', []);
    end

    if colIndex == 0
        ylabel('Log Power Density', 'FontName', 'Arial', 'FontSize', fontSize);
    else
        set(gca, 'YTickLabel', []);
    end

    for i = 1:length(tickMarks)
        xx = tickMarks(i) + [-1 1] * 8E-4 * log(tickMarks(i))^6;
        patch([xx fliplr(xx)], [-20 -20 20 20], '-', 'Facecolor', ones(3, 1) * (0.82 + 0.1 * (i == 5 || i == 1)), 'EdgeColor', 'none');
    end

    if subplotIndex > 1
        plot(timePeriod0, powerSpectrum0, '-', 'color', ones(3, 1) * 0.7, 'LineWidth', 2);
    end

    % Limit data to valid time range
    validIndices = find(tip2 < timeMax);
    numSubVars = size(selectedVariables, 2);

    for varIdx = 1:numSubVars
        currentVar = selectedVariables(1 + subplotIndex, varIdx);

        if currentVar > 1
            dataIndices = find(~isnan(dat(validIndices, currentVar)));
            timeSeriesData = dat(validIndices(dataIndices), currentVar);
            labelStr = legdat{currentVar};
            labelStr = regexprep(labelStr, '_', '');
            labelStr = regexprep(labelStr, '-', '');
            fixedPlot = (currentVar == 35);

            switch currentVar
                case 2
                    labelStr = 'RGR Europe';
                    fixedPlot = 6;
            end

            if currentVar == 2
                plotColor = 'k';
            else
                if subplotIndex <= 1
                    plotColor = lineColors2(subplotIndex * size(selectedVariables, 2) + varIdx, :);
                else
                    plotColor = lineColors(2 * varIdx - 1, :);
                end
            end

            legendLabels{varIdx} = labelStr;

            % FFT Calculation
            samplingFreq = 1E-3 / (tip2(2) - tip2(1));
            numPoints = length(timeSeriesData);
            fftData = fft(timeSeriesData) / numPoints;
            fftData = fftData(1:floor(numPoints / 2) + 1);
            powerDensity = (1 / (samplingFreq * 1)) * abs(fftData).^2;
            powerDensity(2:end-1) = 2 * powerDensity(2:end-1);
            freq = 0:samplingFreq/length(timeSeriesData):samplingFreq/2;
            logPowerDensity = log10(powerDensity);
            period = 1 ./ freq(2:end);
            power = logPowerDensity(2:end);

            interpolatedPower = fixedPlot + interp1(1./period, power, 1./timePeriod);
            smoothedPower = movweighavg(1E5./timePeriod, interpolatedPower, 25, 10);

            % shift down for better comparison
            if mean(smoothedPower) > -1.1
                smoothedPower = smoothedPower - mean(smoothedPower) - 1.2;
            end
            % shift up for better comparison
            if mean(smoothedPower) < -4
                smoothedPower = smoothedPower - mean(smoothedPower) - 1.2;
            end

            if currentVar == 2
                timePeriod0 = timePeriod;
                powerSpectrum0 = smoothedPower;
            end

            fprintf('%d %s\t%1.1f %1.1f\n', varIdx, labelStr, mean(smoothedPower), max(smoothedPower))
            legendHandles(varIdx) = plot(timePeriod, smoothedPower, '-', 'color', plotColor, 'LineWidth', 4 + (currentVar == 2) * 1);
            legendHandles(varIdx).Color(4) = 0.7;
        end
    end
    % clean labels (for publishable plots)
    legendLabels = regexprep(legendLabels, 'PCA', 'PC ');
    legendLabels = regexprep(legendLabels, 'A', ' A');
    legendLabels = regexprep(legendLabels, '\.1', '\.');
    legendLabels = regexprep(legendLabels, 'boombust', 'changes');
    legendLabels = regexprep(legendLabels, 'homogeneity', [newline ' homogeneity']);
    legendLabels = regexprep(legendLabels, '3V', ['3V' newline ' ']);

    % Display legend
    legendHandle = legend(legendHandles, legendLabels);
    set(legendHandle, 'box', 'off', 'fontSize', fontSize - 2, 'location', 'southeast');

    % Add subplot label
    annotation('textbox', [0.084 + colIndex * 1.1 * plotWidth, 0.08 + (rowIndex + 0.79) * 1.05 * plotHeight, .06, .05], ...
        'string', labelLetters(1 + subplotIndex), 'fontsize', 48, 'Fontweight', 'bold', 'EdgeColor', 'none');
end

% Save the plot
outputFile = sprintf('%splots/fft_%d.png', outputDirectory, subplotIndex);
set(gcf, 'PaperPositionMode', 'auto', 'InvertHardCopy', 'off');
print('-dpng', '-r300', outputFile);
