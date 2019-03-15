

function [gsd, d25, d50, d75, d84] = gsdData(dataTable, fileNum, plotFlag)
%
% [gsd, d25, d50, d75, d84] = gsdData(dataTable, fileNum, plotFlag)
%
% gsdData obtains grain size distribution (GSD), calculates some quantiles,
% and optionally plots the GSD.  It requires the 'SurfTableData' or any
% equivalent structured table (see 'headers' cell for column description).
%
% ARGUMENTS:    dataTable = sieving data table mat file
%               fileNum = corresponding line number of 'filepaths' mat file
%               plotFlag = optional logical flag (true produces GSD plot)
%
% RETURNS:      gsd = grain size distribution table (used in plot)
%               d25,...,d84 = percentiles of gsd
%
% NOTES: (1) 'SurfTableData' is structured so that the user can modify this
%             code to return distributions based on multiple criteria (e.g.
%             morphological location, surface/subsurface, or grain size
%             range used in the flume experiments (see headers matrix).
%        (2) 'fileNum' setting to a value of 999 will return GSD for all
%             experimental runs from mass totals taken on surface and
%             subsurface sieving data.
%        (3)  Grain size range is coded for 1.4mm to 63mm sizes only, which
%             corresponds with the coarse grain size limits of the sample
%             image sets (due to resolution limits of the data).
%
% Please see details in the README.md file located on the PATCheS Project
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgContacts, ImgParticles, ImgBedExtend, ImgSurfaces,
% ImgExposure, and ImgEntrainment.

% REFERENCES
% Voepel, H., J. Leyland, R. Hodge, S. Ahmed, and D. Sear (2019),
% Development of a vector-based 3D grain entrainment model with
% application to X-ray computed tomography (XCT)scanned riverbed
% sediment, Earth Surface Processes and Landforms, doi: 10.1002/esp.4608
%
% Copyright (C) 2018  PATCheS Project (http://www.nercpatches.org/)


load(dataTable)

sieveData = data(data(:,30)==fileNum,8:19); % get 63mm to 1.4mm sizes
sieveProb = cumsum(flip(sum(sieveData)))'/sum(sum(sieveData))*100;

%--------INTERPOLATE SIZE FOR 25%TILE--------

% get bounds of sizes
y2 = sizes(sieveProb > 25);
y2 = y2(1);
y1 = sizes(sieveProb <= 25);
y1 = y1(end);

% get bounds of probabilities
x2 = sieveProb(sieveProb > 25);
x2 = x2(1);
x1 = sieveProb(sieveProb <= 25);
x1 = x1(end);

d25 = y1 + (y2 - y1)/(x2 - x1)*(25 - x1);

%--------INTERPOLATE SIZE FOR 50%TILE--------

% get bounds of sizes
y2 = sizes(sieveProb > 50);
y2 = y2(1);
y1 = sizes(sieveProb <= 50);
y1 = y1(end);

% get bounds of probabilities
x2 = sieveProb(sieveProb > 50);
x2 = x2(1);
x1 = sieveProb(sieveProb <= 50);
x1 = x1(end);

d50 = y1 + (y2 - y1)/(x2 - x1)*(50 - x1);

%--------INTERPOLATE SIZE FOR 75%TILE--------

% get bounds of sizes
y2 = sizes(sieveProb > 75);
y2 = y2(1);
y1 = sizes(sieveProb <= 75);
y1 = y1(end);

% get bounds of probabilities
x2 = sieveProb(sieveProb > 75);
x2 = x2(1);
x1 = sieveProb(sieveProb <= 75);
x1 = x1(end);

d75 = y1 + (y2 - y1)/(x2 - x1)*(75 - x1);

%--------INTERPOLATE SIZE FOR 84%TILE--------

% get bounds of sizes
y2 = sizes(sieveProb > 84);
y2 = y2(1);
y1 = sizes(sieveProb <= 84);
y1 = y1(end);

% get bounds of probabilities
x2 = sieveProb(sieveProb > 84);
x2 = x2(1);
x1 = sieveProb(sieveProb <= 84);
x1 = x1(end);

d84 = y1 + (y2 - y1)/(x2 - x1)*(84 - x1);

%-------plotting results-------
if plotFlag

    f = semilogx(sizes,sieveProb);
    f.Marker = 'o';
    f.LineWidth = 3;
    f.Color = 'red';
    xlabel('grain size (mm)')
    ylabel('percent finer (%)')
    ylim([0 100])
    hold on
    plot(sizes,25*ones(size(sizes)),'black')
    plot(sizes,50*ones(size(sizes)),'black')
    plot(sizes,75*ones(size(sizes)),'black')
    plot(sizes,84*ones(size(sizes)),'black')
    set(gca,'FontSize',24)
    hold off

end

gsd = [sizes sieveProb];

end
