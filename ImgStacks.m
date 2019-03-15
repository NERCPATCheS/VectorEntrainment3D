
function ImgStacks(stacksPath, outPath, outFile, idxBottom, resMM)
%
% ImgStacks imports stacked TIFF files of particles and fine-grain matrix.
%
% ImgStacks(stacksPath, outPath, outFile, idxBottom, resMM) is a subroutine
% that requires two input image sets in TIFF stack format: an 8-bit binary
% image of separated sample grains and a 16-bit grayscale of corresponding
% fine-grained sample matrix obrained during image segmentation (see Wiki).
%
% ImgStacks requires the MATLAB Image Processing Toolbox, the Statistics
% and Machine Learning Toolbox, and has the following subroutine arguments:
%
%   stacksPath = path from current working directory to TIFF stack folders
%   outPath = path from current working directory to write all output images
%   outFile = name of MAT file to store output image arrays and metrics
%   idxBottom = index of 2D layer that contains bottom of the scanned sample
%   resMM = voxel side length resolution for the scanned sample image (mm)
%
% ImgStacks writes the following TIFF virtual stacks as binary (bw) and
% labelled (label) images in the 'outPath' sample subfolder:
%
%   bwMatrix = binary image of fine-grained matrix threshold at mean grayscale
%   bwParticles = binary image of whole particles (excluding cropped grains)
%   bwParticlesFull = binary image of all particles (same as input stack)
%   labelParticles = same as 'labelParticlesFull' with cropped grains removed
%   labelParticlesFull = labelled grains used as ID (ordered top to bottom)
%   labelParticlesSurf = surface voxels for all grains in 'labelParticlesFull'
%
% ImgStacks saves these same binary and labelled 3D images as arrays along
% with their corresponding connected component (cc) structure arrays to the
% 'outFile' MAT file saved in the current working directory.  The structure
% array 'dataParticles' is created and saved to 'outFile', and contains the
% following metrics for each particle in the 'labelParticlesFull' image:
%
%   ParticleID = unique grain identifier from the 'labelParticleFull' array
%   Centroid = 3D coordinates of the grain centre of mass (uniform density)
%   VolumeMM3 = grain volume (mm3) found from image voxel resolution (in mm)
%   BoundingBox = grain bounding box used to find local mean bed elevation
%   CroppedFlag = logical (whole stone = 0, stone cropped off at bottom = 1)
%
% Please see details in the README.md file located on the PATCheS Project
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgContacts, ImgParticles, ImgBedExtend, ImgSurfaces, ImgExposure,
% and ImgEntrainment.

% REFERENCES
% Voepel, H., J. Leyland, R. Hodge, S. Ahmed, and D. Sear (2019),
% Development of a vector-based 3D grain entrainment model with
% application to X-ray computed tomography (XCT)scanned riverbed
% sediment, Earth Surface Processes and Landforms, doi: 10.1002/esp.4608
%
% Copyright (C) 2018  PATCheS Project (http://www.nercpatches.org/)

%---------CHECKING REQUIREMENTS BEFORE RUN------------

% Check user has required toolbox(s) installed installed
tbCheck = license('test', 'Image_Toolbox') && license('test', 'Statistics_Toolbox');
if ~tbCheck
	% User does not have the toolbox installed.
	error('Requires Image Processing Toolbox & Statistics and Machine Learning Toolbox.')
end


%---------IMPORT BINARY PARTICLE IMAGE STACK----------


% import tiff filenames for stacks
pathParticlesFull = char(strcat(pwd,stacksPath,'ParticlesFull/'));
tifParticlesFull = dir([pathParticlesFull '*.tif']);

% get information on tiffs
fileNameParticlesFull = [pathParticlesFull tifParticlesFull(1).name];
infoParticlesFull = imfinfo(fileNameParticlesFull);
tifCountFull = length(tifParticlesFull); % use for image tiff stacks

% get image dimensions
imgWidthFull = infoParticlesFull(1).Width;
imgHeightFull = infoParticlesFull(1).Height;

% preallocating logic space for image
bwParticlesFull = false(imgHeightFull,imgWidthFull,tifCountFull);

% retrieve image slices from tiff stack
for k = 1:tifCountFull
    fileNameParticlesFull = [pathParticlesFull tifParticlesFull(k).name];
    bwParticlesFull(:,:,k) = imread(fileNameParticlesFull);
end



%---------IMPORT GRAYSCALE MATRIX IMAGE STACK----------

% import tiff filenames for stacks
pathMatrix = char(strcat(pwd,stacksPath,'Matrix/'));
tifMatrix = dir([pathMatrix '*.tif']);

% get information on tiffs
fileNameMatrix = [pathMatrix tifMatrix(1).name];
infoMatrix = imfinfo(fileNameMatrix);
tifCount = length(tifMatrix); % use for image tiff stacks

% get image dimensions
imgWidth = infoMatrix(1).Width;
imgHeight = infoMatrix(1).Height;

% preallocating unsigned 16bit space for image
gsMatrix = uint16(zeros(imgHeight,imgWidth,tifCount));

% retrieve image slices from tiff stack
for k = 1:tifCount
    fileNameMatrix = [pathMatrix tifMatrix(k).name];
    gsMatrix(:,:,k) = imread(fileNameMatrix);
end

% performing grayscale open to enhance thresholding
se = strel3d(3);
gsMatrix = imopen(gsMatrix,se);

% getting mean of 16bit grayscale image stack without zero values
tbl = tabulate(gsMatrix(:));
tbl(1,:) = []; % removing zero counts
one = ones(length(tbl),1); % ones vector to find mean
mu = tbl(:,1)'*tbl(:,2)/(one'*tbl(:,2)); % finding mean

% generating binary matrix thresholding at mean value
bwMatrix = false(imgHeight,imgWidth,tifCount);
for k = 1:tifCount
    fileNameMatrix = [pathMatrix tifMatrix(k).name];
    bwMatrix(:,:,k) = imbinarize(imread(fileNameMatrix), mu/2^16);
end

% trimming any matrix overlap with stones
bwPFullComp = imcomplement(bwParticlesFull);
bwMatrix = immultiply(bwPFullComp,bwMatrix);



%---------MAKING LABELLED PARTICLES AND SURFACES----------

% getting connectivity of binary amd making labelled particles
connParticles = 6;
ccParticlesFull = bwconncomp(bwParticlesFull,connParticles);
labelParticlesFull = uint16(labelmatrix(ccParticlesFull));

% make label particle surfaces (hollow shells of stones)
bwPFullErod = imerode(bwParticlesFull,se);
bwPFullErodComp = imcomplement(bwPFullErod);
labelParticlesSurf = immultiply(bwPFullErodComp,labelParticlesFull);



%---------REMOVING CROPPED LABELLED PARTICLES----------

% remove cropped stones from labelled particles
bottomParticles = unique(labelParticlesFull(:,:,idxBottom:end));
includeStoneIDs = setdiff(1:ccParticlesFull.NumObjects,bottomParticles);

% preallocating binary to built particle set of whole stones
bwParticles = false(size(bwParticlesFull));

% construct image that excludes cropped stones
for i = includeStoneIDs

    bwParticles(ccParticlesFull.PixelIdxList{i}) = true;

end

labelParticles = immultiply(labelParticlesFull,bwParticles);




%---------GETTING GEOMETRY FOR EACH PARTICLE----------

% preallocate storage for Particle-to-Contacts connection metrics
dataParticles = repmat(struct(...
    'ParticleID',0,'Centroid',-1,'VolumeMM3',0,'BoundingBox',0,...
    'CroppedFlag',0),ccParticlesFull.NumObjects,1);

% get centroids, volumes and bounding boxes of particles
% note: Area is by pixel count in 2D; so, a 3D voxel count is volume
part = regionprops(ccParticlesFull,'Centroid','BoundingBox','Area');

% storing information for particles
for k = 1:ccParticlesFull.NumObjects

    dataParticles(k).ParticleID = k;
    dataParticles(k).Centroid = part(k).Centroid; %([2 1 3]); % convert
    dataParticles(k).VolumeMM3 = part(k).Area*resMM^3;
    dataParticles(k).BoundingBox = part(k).BoundingBox;

    % flag cropped stones
    if isempty(bottomParticles(bottomParticles==k))
        dataParticles(k).CroppedFlag = false;
    else
        dataParticles(k).CroppedFlag = true;
    end

end



%---------SAVING TO MAT AND WRITING IMAGES TO OUT FOLDERS----------


% save arrays for analysis in other scripts
save(outFile,'ccParticlesFull','bwParticles','bwParticlesFull','bwMatrix',...
    'labelParticles','labelParticlesFull','labelParticlesSurf','dataParticles')


% Write BW and Label Images to 3D Multi-Layered TIFF Image files
fPart = char(strcat(pwd,outPath,'labelParticles.tif'));
fPartFull = char(strcat(pwd,outPath,'labelParticlesFull.tif'));
fPartSurf = char(strcat(pwd,outPath,'labelParticlesSurf.tif'));
bwPart = char(strcat(pwd,outPath,'bwParticles.tif'));
bwPartFull = char(strcat(pwd,outPath,'bwParticlesFull.tif'));
bwMat = char(strcat(pwd,outPath,'bwMatrix.tif'));

imwrite(uint16(labelParticles(:,:,1)),fPart,'tif','Compression','none')
for k = 2:tifCountFull
    imwrite(uint16(labelParticles(:,:,k)),fPart,'tif',...
        'Compression','none','WriteMode','append')
end

imwrite(uint16(labelParticlesFull(:,:,1)),fPartFull,'tif','Compression','none')
for k = 2:tifCountFull
    imwrite(uint16(labelParticlesFull(:,:,k)),fPartFull,'tif',...
        'Compression','none','WriteMode','append')
end

imwrite(uint16(labelParticlesSurf(:,:,1)),fPartSurf,'tif','Compression','none')
for k = 2:tifCountFull
    imwrite(uint16(labelParticlesSurf(:,:,k)),fPartSurf,'tif',...
        'Compression','none','WriteMode','append')
end

imwrite(uint16(bwParticles(:,:,1)),bwPart,'tif','Compression','none')
for k = 2:tifCountFull
    imwrite(uint16(bwParticles(:,:,k)),bwPart,'tif',...
        'Compression','none','WriteMode','append')
end

imwrite(uint16(bwParticlesFull(:,:,1)),bwPartFull,'tif','Compression','none')
for k = 2:tifCountFull
    imwrite(uint16(bwParticlesFull(:,:,k)),bwPartFull,'tif',...
        'Compression','none','WriteMode','append')
end

imwrite(uint16(bwMatrix(:,:,1)),bwMat,'tif','Compression','none')
for k = 2:tifCount
    imwrite(uint16(bwMatrix(:,:,k)),bwMat,'tif',...
        'Compression','none','WriteMode','append')
end

end % end imgStacks function

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
