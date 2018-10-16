
function ImgBedExtend(outFile, outPath, radTrim)
%
% ImgBedExtend extends image sample edge to mimic natural bed of channel.
%
% ImgBedExtend(outFile, outPath, radTrim) is a subroutine that replicates
% binary copies of the scanned sample image around the edge of the original
% sample image. This effectively extends the bed around the original sample
% image so that hiding effects of upstream stones may be realised.  Where 
% any stone exists, the binary logical is true. It is then converted to a
% label matrix where all stones in the extended bed are assigned a common
% ParticleID = n + 1 for n total stones in the sample (see README).
%
% ImgBedExtend requires the MATLAB Image Processing Toolbox, the Mapping
% Toolbox, and has the following subroutine arguments:
% 
%   outFile = name of MAT file to store a 3D label array of the extended bed
%   outPath = path from current working directory to write an output image
%   radTrim = a voxel distance to overlap original sample before trimmming
% 
% ImgBedExtend writes the final 'labelParticlesBed' to a 3D TIFF image in
% the 'outPath' folder and saves the 3D image array in the 'outFile'.
%
% Please see details in the README.md file located on the PATCheS Project 
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgContacts, ImgParticles, ImgSurfaces, ImgExposure,
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
tbCheck = license('test', 'Image_Toolbox') && license('test', 'Map_Toolbox');
if ~tbCheck
	% User does not have the toolbox installed.
	error('Requires Image Processing Toolbox & Mapping Toolbox.')
end


%---------FIND XYZ-CENTER OF PARTICLES AND RADII FOR COPIES----------

% load data, get image size and particle count
load(outFile)
[x, y, z] = size(labelParticlesFull);
n = length(dataParticles);   

% get coordinate bounds from BoundingBox field
bb = extractfield(dataParticles,'BoundingBox');
bb = reshape(bb,6,length(bb)/6)';

% find XYZ means of smallest/largest 10 elements
bbMin = sort(bb(:,1:3));
bbMax = sort(bb(:,1:3)+bb(:,4:6),'descend');
minXYZ = floor(mean(bbMin(1:10,:)));
maxXYZ = ceil(mean(bbMax(1:10,:)));

% calculate XYZ coordinate centroids
centX = floor(mean([minXYZ(1) maxXYZ(1)]));
centY = floor(mean([minXYZ(2) maxXYZ(2)]));
centZ = floor(mean([minXYZ(3) maxXYZ(3)]));
centXYZ = [centX centY centZ];

% make trimmed XY radius to overlap labelled stones 
R = max(maxXYZ(1:2) - minXYZ(1:2)) - radTrim;

% make spherical radius to trim overall image
dXYZ = maxXYZ - centXYZ;
Rtrm = sqrt(dXYZ*dXYZ') + 10;



% ---------GENERATE BINARY BED, ADD MATRIX AND THEN TRIMMED------------

% making binary copies of particles to surround labelled particles
rot = 0:30:360; % rotation angles for repeated copies
xyRot = [R*cos(rot(1)*pi/180) R*sin(rot(1)*pi/180)];
bwBed = imtranslate(bwParticlesFull,round(xyRot));
for k = 2:length(rot)
    xyRot = [R*cos(rot(k)*pi/180) R*sin(rot(k)*pi/180)];
    bwBed = bwBed | imtranslate(bwParticlesFull,round(xyRot));
end

% adding matrix to bed and then trim to particles
bwBed = bwBed | bwMatrix;
bwTrim = imdilate(bwParticlesFull,strel3d(3));
bwTrim = imcomplement(bwTrim);
bwBed = bwBed & bwTrim;

% make binary sphere to trim new binary image boundary
[X, Y, Z] = meshgrid(1:y, 1:x, 1:z); % coords must swap due to rotation
trmSphere = sqrt((X-centX).^2+(Y-centY).^2+(Z-centZ).^2) <= Rtrm;
bwBed = bwBed & trmSphere;

% assign new ID value to binary wall and matrix then add to particles
wallID = n + 1; % incrementing ID value
lblBed = uint16(wallID*bwBed); % tag wall with identification
labelParticlesBed = imadd(lblBed, labelParticlesFull); % adding wall



% ---------SAVING PARTICLES FULL AND WRITING TO TIFF FILE--------------

% write combined image to file for inspection
lblFile = char(strcat(pwd,outPath,'labelParticlesBed.tif')); %make filename

imwrite(labelParticlesBed(:,:,1),lblFile,'tif','Compression','none')
for k = 2:z
    imwrite(labelParticlesBed(:,:,k),lblFile,'tif',...
        'Compression','none','WriteMode','append')
end

% save label particles image containing random wall
save(outFile,'labelParticlesBed','-append')

end % end ImgBedExtend function


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

