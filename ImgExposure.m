
function ImgExposure(outFile, viewAngles, viewImages, surfaces, resMM)
%
% ImgExposure calculates exposure, areas, and a frontal image for each grain.
%
% ImgExposure(outFile, viewAngles, viewImages, surfaces, resMM) uses two
% cell vectors of images returned by the 'ImgSurface' function--'viewImages'
% and 'surfaces'--to calculate an exposure ratio vector at 'viewAngles', a
% corresponding exposure factor, overhead and frontal cross-sectional grain 
% areas, and a frontal view binary image of a single grain for each stone.
% 
% ImgExposure requires the MATLAB Parallel Computing Toolbox, the MATLAB
% Image Processing Toolbox, the MATLAB Statistics and Machine Learning 
% Toolbox and calls on the following subroutine arguments:
% 
%   outFile = name of MAT file to store the updated 'dataParticles' array
%   viewAngles = set of bed viewing angles where zero is the overhead view
%   viewImages = cell vector of 3D labelled images rotated at 'viewAngles'
%   surfaces = cell vector of 2D projected images for each 'viewImages'
%   resMM = voxel side length resolution for the scanned sample image (mm)
% 
% ImgExposure updates and saves the 'dataParticles' structure array with
% the following metrics for each stone in the sample:
%
%   ExposureRatio = exposed-to-total 2D grain area ratio at each 'viewAngles'
%   ExposureFactor = calculated from each 'ExposureRatio' at each 'viewAngles'
%   FrontalAreaMM2 = grain cross-sectional area perpendicular to streamflow (mm)
%   OverheadAreaMM2 = grain cross-sectional area parellel to streamflow (mm)
%   FrontalImage = binary 2D projected image of grain viewed from upstream
%
% ImgExposure saves the 'dataParticles' and 'viewAngles' to the 'outFile'.
%
% Please see details in the README.md file located on the PATCheS Project 
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgContacts, ImgParticles, ImgBedExtend, ImgSurfaces,
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
tbCheck = license('test', 'Image_Toolbox') ...
    && license('test', 'Distrib_Computing_Toolbox') ...
    && license('test', 'Statistics_Toolbox');
if ~tbCheck
	% User does not have the toolbox installed.
	error('Requires Image Processing, Parellel Computing & Statistics and Machine Learning Toolboxes.')
end

%---------USE SURFACES TO CALCULATE EXPOSURE----------

tic

load(outFile)

% tabulate surface labels, calculate exposure ratios, and construct table
fprintf('Calculating exposure ratios for %s...\n',outFile)

% getting number of stones and viewing angles
n = ccParticlesFull.NumObjects; % getting stone count
m = length(viewAngles); % get total number of angles used
fprintf('Total stone count = %s\n',num2str(n))

% getting indices for overhead and downstream view angles
k0 = find(viewAngles == 0); % index for overhead (0deg) viewing image
k90 = find(viewAngles == 90); % index for downstream (90deg) viewing image

% preallocating for exposure ratios, 2D planar areas, and downstream image
surfaceTable = zeros(n,m); % tabulated exposure ratio vectors for each stone
frontalArea = zeros(n,1); % vector to store frontal areas (drag force)
overheadArea = zeros(n,1); % vector to store overhead areas (lift force)
frontalImage = cell(n,1); % binary image of downstream view of single stone



for k = 1:m % loop over surfaces (excluding 270 degree surface)

    % tabulate surface label values
    tbl = tabulate(surfaces{k}(:));
    tbl(1,:) = []; % dropping label = 0 (background)
    tbl(end,:) = []; % dropping extended bed (ParticleID = n + 1)

    % get linear pixel index list of current label array
    labels = regionprops(viewImages{k},'PixelIdxList');

    % assigning kth view image for parfor processing
    viewImagesK = viewImages{k};

    % get exposure ratios
    parfor i = 1:n % loop over stones

        % construct image for single stone, ParticleID = i
        particle = false(size(viewImagesK));
        particle(labels(i).PixelIdxList) = true;

        % union all slices to get area of stone viewed at kth angle
        imSlice = particle(:,:,1);
        for j = 2:size(particle,3)
            imSlice = imSlice | particle(:,:,j);    
        end

        % storing downstream and overhead viewed areas (mm2)
        if k == k90 % downstream area (for drag force)
            frontalArea(i) = sum(imSlice(:))*resMM^2; % in mm2
            frontalImage{i} = imSlice; % storing downstream area image
        elseif k == k0 % overhead area (for lift force)
            overheadArea(i) = sum(imSlice(:))*resMM^2; % in mm2
        end

        % tabulate exposure ratio and frontal area
        if isempty(tbl(tbl(:,1)==i,2))
            surfaceTable(i,k) = 0; % if ith stone not in table
        else % calculate exposure ratio vector
            surfaceTable(i,k) = tbl(tbl(:,1)==i,2)/sum(imSlice(:));
        end

    end % end loop over stones

    fprintf('Finished with exposure ratios for angle = %s\n',num2str(viewAngles(k)))

end % end loop over surfaces

% killing parallel workers
poolobj = gcp('nocreate');
delete(poolobj);

% calculating exposure factors for each particle
fprintf('Calculating exposure factors for %s\n',outFile)

% preallocating storage space
dataParticles(n).ExposureRatio = 0;
dataParticles(n).ExposureFactor = 0;
dataParticles(n).FrontalAreaMM2 = 0;
dataParticles(n).OverheadAreaMM2 = 0;
dataParticles(n).FrontalImage = 0;

% normalizing sine of angles such that weights sum to unity
normSin = sin(viewAngles*pi/180).^2/sum(sin(viewAngles*pi/180).^2);

for i = 1:n
    dataParticles(i).ExposureRatio = surfaceTable(i,:);
    dataParticles(i).ExposureFactor = normSin*surfaceTable(i,:)';
    dataParticles(i).FrontalAreaMM2 = frontalArea(i);
    dataParticles(i).OverheadAreaMM2 = overheadArea(i);
    dataParticles(i).FrontalImage = frontalImage{i};
end

fprintf('Updating dataParticles mat file %s\n',outFile)
save(outFile,'dataParticles','viewAngles','-append')

toc

end % end function


