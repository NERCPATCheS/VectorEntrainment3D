
function [viewImages, surfaces] = ImgSurfaces(outFile, outPath, viewAngles)
%
% ImgSurfaces generates a set of tilted 3D images and 2D projected surfaces.
%
% [viewImages, surfaces] = ImgSurfaces(outFile, outPath, viewAngles) uses a
% set of user defined viewing angles from 0 to 90 degrees (see NOTE below)
% to generate tilted 3D images and their corresponding 2D projected surfaces
% for use in subsequent subroutine calculations.
%
% NOTE: User must include both 0 and 90 degree angles as these are used in
% the subsequent 'ImgExposure' subroutine run to calculate and store both
% parallel and perpendicular 2D cross-sectional grain areas and to generate
% a downstream image of each grain used in drag force calculations.
%
% ImgSurfaces requires the MATLAB Image Processing Toolbox, and has the
% following function arguments:
%
%   outFile = name of MAT file to store metrics for each particle
%   outPath = path from current working directory to write tilted bed images
%   viewAngles = vector of bed viewing angles where zero is the overhead view
%
% ImgSurfaces creates two folders in the 'outPath' folder, 'Rotated/' and
% 'Surface/', to write 3D rotated images and 2D rotated label surfaces as
% TIFF images under path and filenames
%
%   <outPath>Rotated/labelParticlesBedXXX.tif, and
%   <outPath>Surface/labelSurfaceXXX.tif,
%
% where XXX is the user defined viewing angle (bewtween 0 and 90 degrees).
%
% ImgSurfaces returns two cell vectors of the same length as 'viewAngles'
% where each cell element contains a rotated image as follows:
%
%   viewImages = 3D rotated images of 'labelPartilesBed' for each angle
%   surfaces = 2D projected images of each 'viewImages' for each angle
%
% IMPORTANT: The 'ImgExposure' subroutine uses both returned vectors in
% further calculations; therefore, 'ImgSurfaces' and 'ImgExposure' should
% be run sequentially for each sample so as not to overwrite these vectors.
%
% Please see details in the README.md file located on the PATCheS Project
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgContacts, ImgParticles, ImgBedExtend, ImgExposure,
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
tbCheck = license('test', 'Image_Toolbox');
if ~tbCheck
	% User does not have the toolbox installed.
	error('Requires Image Processing Toolbox.')
end

%---------GENERATING ROTATED IMAGE SURFACES-----------

load(outFile)

% STEP 1: Make angled views from labelParticlesWall image

% generate 3D affine rotational matrices for all viewing angles
tform = RotationMatrix(viewAngles);
n = length(tform);

% make filepaths with degree angles for saving 3D rotated images
childFolder = 'Rotated/';
fileRoot = 'labelParticlesBed';
rotFilename = MakeFileNames(outPath, fileRoot, childFolder, viewAngles);

% making rotated images of labelled particles with wall
viewImages = cell(size(tform)); % preallocating for images
for i = 1:n

    % rotate images through all angles
    viewImages{i} = imwarp(labelParticlesBed,tform{i},'nearest');
    tiffCount = size(viewImages{i},3);

    % write rotated 3D image to file
    imwrite(uint16(viewImages{i}(:,:,1)),rotFilename{i},'tif','Compression','none')
    for k = 2:tiffCount
        imwrite(uint16(viewImages{i}(:,:,k)),rotFilename{i},'tif',...
            'Compression','none','WriteMode','append')
    end

end

% STEP 2: Get 2D surfaces for each viewAngle plus 270 degree image
disp(strcat('Generating surface label images for', sprintf(' %s',outFile)))

% preallocating for surfaces
surfaces = cell(size(viewImages));

% make filepaths with degree angles for saving surface images
childFolder = 'Surface/';
fileRoot = 'labelSurface';
surfFilename = MakeFileNames(outPath, fileRoot, childFolder, viewAngles);

% make surface images
for k = 1:n % loop over image views

    % getting size of array
    [x, y, z] = size(viewImages{k});

    surfaces{k} = uint16(zeros(x,y));
    for i = 1:x
        for j = 1:y

            % get z-vector of labelParticles at location (x,y)
            zlabPart = reshape(viewImages{k}(i,j,:),z,1);

            % get first nonzero element (which is a labelled surface)
            zlabPart = zlabPart(zlabPart > 0);

            if isempty(zlabPart)
                surfaces{k}(i,j) = 0;
            else
                surfaces{k}(i,j) = zlabPart(1);
            end

        end
    end

    % write 2D surface images to tiff file
    imwrite(uint16(surfaces{k}),surfFilename{k},'tif','Compression','none')

end % end image view loop

end % end function


%=======================INTERNAL FUNCTIONS=============================

% makes filenames used to write 3D rotations and their 2D surfaces
function [filename] = MakeFileNames(outPath, fileRoot, childFolder, angles)

    % make filenames to write surface tiffs
    filename = cell(size(angles));
    outFolder = char(strcat(pwd,outPath,childFolder));

    % make folder if it doesn't exist
    if not(7 == exist(outFolder,'dir'))
        mkdir(outFolder)
    end

    for i = 1:length(filename)

        % pad filenames with zeros
        if angles(i) < 10
            filename{i} = strcat(outFolder, fileRoot, '00', num2str(angles(i)), '.tif');
        elseif angles(i) < 100
            filename{i} = strcat(outFolder, fileRoot, '0', num2str(angles(i)), '.tif');
        else
            filename{i} = strcat(outFolder, fileRoot, num2str(angles(i)), '.tif');
        end

    end

end


% makes affine 3D rotation matrices for all angles (out is cell)
function [tform] = RotationMatrix(angles)

    theta = pi/180*angles;
    n = length(theta);
    c = cos(theta);
    s = sin(theta);
    tform = cell(size(theta));
    for i = 1:n

        ty = [c(i) 0 s(i) 0
                0  1   0  0
             -s(i) 0 c(i) 0
                0  0   0  1];
        tform{i} = affine3d(ty);

    end

end
