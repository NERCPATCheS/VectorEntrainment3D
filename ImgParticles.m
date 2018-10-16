
function ImgParticles(outFile, resMM)
%
% ImgParticles calculates granular metrics for each particle in the sample.
%
% ImgParticles(outFile, resMM) is a subroutine that calculates principal
% axis lengths and their corresponding spatial orientation, the maximum 
% elevation of the particle and its coordinates, and the fraction of total
% granular surface area in contact with the fine-grained matrix. 
%
% ImgParticles requires the MATLAB Image Processing Toolbox, the Statistics 
% and Machine Learning Toolbox, and has the following subroutine arguments:
% 
%   outFile = name of MAT file to store granular metrics for each stone
%   resMM = voxel side length resolution for the scanned sample image (mm)
% 
% ImgParticles appends an update to the structure array 'dataParticles' that
% is saved to the 'outFile' MAT file and contains the following additional 
% metrics for each particle:
%
%   MajorAxisMM = the major principal axis of the particle (mm)
%   MedianAxisMM = the median principal axis of the particle (mm)
%   MinorAxisMM = the minor principal axis of the particle (mm)
%   AxisVectors = a 3x3 matrix that describes particle spatial orientation
%   MaxElevMM = the maximum elevation of the particle in the sample (mm)
%   MaxElevCoord = the 3D coordinates of MaxElevMM as voxel element
%   MatrixContactAreaRatio = matrix fines contact area-t0-total area ratio
%
% Metrics for cropped off stones are empty sets for axis lengths and their 
% corresponding spatial orientation.  AxisVectors are an orthonormal set of
% 3D unit vectors that define the spatial orientation of the principal axes.  
% MaxElevMM has lower values for upper most grains and higher values for
% lower most grains since the z-axis is in the direction of gravity.
%
% Please see details in the README.md file located on the PATCheS Project 
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgContacts, ImgBedExtend, ImgSurfaces, ImgExposure, 
% and ImgEntrainment.

% REFERENCES 
% Voepel, H., J. Leyland, R. Hodge, S. Ahmed, and D. Sear (submitted), 
% Development of a vector-based 3D grain entrainment model with 
% application to X-ray computed tomography (XCT)scanned riverbed
% sediment, Earth Surface Processes and Landforms (?????)
% 
% Copyright (C) 2018  PATCheS Project (http://www.nercpatches.org/)

%---------CHECKING REQUIREMENTS BEFORE RUN------------

% Check user has required toolbox(s) installed installed
tbCheck = license('test', 'Image_Toolbox') && license('test', 'Statistics_Toolbox');
if ~tbCheck
	% User does not have the toolbox installed.
	error('Requires Image Processing Toolbox & Statistics and Machine Learning Toolbox.')
end

%---------PROCESSING MATRIX FOR SURFACE CONTACT----------

tic % checking elapse time

load(outFile)


% open binary matrix to trim flex
se = strel3d(3);
bwMatrixOpen = imopen(bwMatrix,se);

% make table of particle surface areas
tblPartArea = tabulate(labelParticlesSurf(:));
tblPartArea(1,:) = []; % removing zeros

% use dilated binary matrix to grab nearby particle areas
bwMatrixDila = imdilate(bwMatrixOpen,se);
labelPFullCohe = immultiply(bwMatrixDila,labelParticlesSurf);

% make table of particle surface areas in contact with matrix
tblPartCohe = tabulate(labelPFullCohe(:));
tblPartCohe(1,:) = []; % removing zeros


%----------GETTING PRINCIPAL AXES OF PARTICLES-----------

% preallocate space for axes
% NOTE: Please use whole stone counts
dataCount = ccParticlesFull.NumObjects;
dataParticles(dataCount).MajorAxisMM = 0;
dataParticles(dataCount).MedianAxisMM = 0;
dataParticles(dataCount).MinorAxisMM = 0;
dataParticles(dataCount).AxisVectors = 0;
dataParticles(dataCount).MaxElevMM = 0;
dataParticles(dataCount).MaxElevCoord = 0;
dataParticles(dataCount).MatrixContactAreaRatio = 0;

% % for checking stone rotations (comment when finished)
% stoneCheck = zeros(dataCount,3);

for k = 1:dataCount

    % get pixel list and preallocate xyz storage
    elm = ccParticlesFull.PixelIdxList{k};

    % make image of single stone and get size
    stone = false(size(bwParticlesFull));
    stone(elm) = true;
    [x, y, z] = ind2sub(size(stone), find(stone == 1));
    xyzStone = [x, y, z];

    % get pca of xyz
    [xyzEigen, xyzRot, ~] = pca(xyzStone);

    % get principal axes and original xyz lengths in milimeters
    principalAxes = (max(xyzRot) - min(xyzRot))*resMM;
    
%     % for checking stone rotations (comment when finished)
%     stoneCheck(k,:) = (max(xyzStone) - min(xyzStone))*resMM;

    % get maximum protrusion height and its coordinate indices
    dataParticles(k).MaxElevMM = min(xyzStone(:,3))*resMM;
    maxProtrusions = xyzStone(xyzStone(:,3)==min(xyzStone(:,3)),:);

    % check whether a mean of protrusions is needed (note: swapped XY)
    if isvector(maxProtrusions)
        maxElevCoord = round(maxProtrusions);
    else
        maxElevCoord = round(mean(maxProtrusions));
    end

    % must swap XY coordinates using [2 1 3] due to rotation
    dataParticles(k).MaxElevCoord = maxElevCoord([2 1 3]);

    % store prinicpal axes lengths
    if length(principalAxes) < 3 || dataParticles(k).CroppedFlag
        dataParticles(k).MajorAxisMM = [];
        dataParticles(k).MedianAxisMM = [];
        dataParticles(k).MinorAxisMM = [];
        dataParticles(k).AxisVectors = [];
    else
        dataParticles(k).MajorAxisMM = principalAxes(1);
        dataParticles(k).MedianAxisMM = principalAxes(2);
        dataParticles(k).MinorAxisMM = principalAxes(3);
        dataParticles(k).AxisVectors = xyzEigen;
    end

    % storing matrix contact metrics
    try
        totalArea = tblPartArea(tblPartArea(:,1)==k,2);
        partArea = tblPartCohe(tblPartCohe(:,1)==k,2);
        dataParticles(k).MatrixContactAreaRatio = partArea/totalArea;
    catch ME
        dataParticles(k).MatrixContactAreaRatio = 0;
    end

    if isempty(dataParticles(k).MatrixContactAreaRatio)
        dataParticles(k).MatrixContactAreaRatio = 0;
    end

    fprintf('...finished with ParticleID = %s, Coord = [%s]\n',...
        num2str(k), num2str(maxElevCoord([2 1 3])))


end


% add updated 'dataParticles' to 'Run3TailArrays.mat' file
save(outFile,'dataParticles', '-append')

toc

end


%====================FUNCTION BLOCK========================

function [se] = strel3d(sz)

% modify strel for 3D spherical SE
sw = (sz - 1)/2; 
ses2 = ceil(sz/2); 
[y,x,z] = meshgrid(-sw:sw, -sw:sw, -sw:sw); 
m = sqrt(x.^2 + y.^2 + z.^2); 
b = (m <= m(ses2, ses2, sz)); 
se = strel('arbitrary', b);

end
