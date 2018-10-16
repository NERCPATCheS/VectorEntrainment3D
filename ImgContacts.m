
function ImgContacts(outFile,seDiam,sampleRate,outGroup)
%
% ImgContacts obtains contact points between object stone and other grains.
%
% ImgContacts(outFile,seDiam,sampleRate,outGroup) is a subroutine that finds
% all contact points between each object stone and its neighbouring grains, 
% and appends results into the 'dataParticles' structure array.
%
% ImgContacts requires the MATLAB Parallel Computing Toolbox, the MATLAB
% Image Processing Toolbox, the MATLAB Statistics and Machine Learning 
% Toolbox and calls on the following subroutine arguments:
% 
%   outFile = name of MAT file to store the Nx8 contact stone metrics
%   seDiam = voxel distance around object stone to 'grab' local contact stones
%   sampleRate = fraction of object stone surface coordinates to resample
%   outGroup = object stone ID to write local group surfaces to TIFF file
% 
% ImgContacts assigns an Nx8 matrix called 'Contacts' as a metric to each
% object stone in 'dataParticles', where N is the number of neighbouring
% particles in proximal contact with the object stone, and the columns are:
%
%   Column 1 is ParticleID of the contact stone to the current object stone
%   Column 2 is squared voxel distance between the object and contact stone
%   Columns 3-5 are the 3D coordinates of the proximal contact point
%   Columns 6-8 are the 3D coordiantes of the contact stone centroid
%
% The structure array 'dataParticles' is updated in the 'outFile' MAT file.
%
% WARNING: Too high of a 'sampleRate' value may result in computer freeze
% due to lack of resources.  Lowering the this value may be necessary.
%
% Please see details in the README.md file located on the PATCheS Project 
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgParticles, ImgBedExtend, ImgSurfaces, ImgExposure,
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


%-------------SEARCH FOR CONTACT POINTS---------------

tic % checking elapse time

load(outFile)

% preallocating space within dataParticles
dataCount = ccParticlesFull.NumObjects;
dataParticles(dataCount).Contacts = 0;

for k = 1:dataCount

    % dat = dataParticles(k)

    % make group of labelled surfaces nearby stone k
    bwStone = labelParticlesFull == k;
    bwStoneDial = imdilate(bwStone,strel3d(seDiam)); % dilate as needed to grab group
    lblStones = immultiply(uint16(bwStoneDial),labelParticlesSurf); % make group

    % write selected labelled surface group to file (default, outGroup = 0 writes nothing)
    if k == outGroup
        imwrite(lblStones(:,:,1),'lblStones15.tif','tif','Compression','none')
        for slc = 2:size(lblStones,3)
            imwrite(lblStones(:,:,slc),'lblStones15.tif','tif',...
                'Compression','none','WriteMode','append')
        end
    end

    % get voxel counts for each surface in group
    tab = tabulate(lblStones(:));
    tab(1,:) = []; % removing background (label = 0) from list

    % get coordinates for each surface in group
    cnt = size(tab,1);
    coord = cell(cnt,1);
    parfor j = 1:cnt
        [x, y, z] = ind2sub(size(lblStones), find(lblStones == tab(j,1)));
        coord{j} = [y, x, z]; % convert from [row col layer] to [x y z]
    end


    %------------get coords and their voxel count for stone k-------------

    % reassign kth object stone coordinates and remove it from table
    idx = find(tab(:,1)==k); 
    x = coord{idx}; % reassign kth stone coordinates as x
    coord(idx,:) = []; % removing kth coordinates from coord
    tab(idx,:) = []; % removing kth entry from table
    n = size(x,1); % voxel count for object stone

    % resample kth stone with fewer coordinates if too large
    if n >= max(tab(:,2)) % get at least sampleRate*100% resample 
        sampleSize = max(round(sampleRate*n),max(tab(:,2)));
        idx = sort(datasample(1:n,sampleSize,'Replace',false));
        x = x(idx,:); % reduce kth coordinate vector
        n = size(x,1); % get new voxel count
    end


    % preallocating for contact list
    cnt = size(tab,1);
    contact = zeros(cnt,8);


    % get coords for nearby stone
    for j = 1:cnt

        % get coordinates list for jth contact point and its length
        yj = coord{j};
        m = size(yj,1);

        % preallocate for minimum squared distances for each yj coord
        dist2min = zeros(m,4);

        % get distances between kth stone and jth stone then pick smallest
        parfor t = 1:m % looping over each voxel
            y = repelem(yj(t,:),n,1); % vector of voxel t length(x)
            dx = x - y;
            dist2 = diag(dx*dx');
            midCoord = (x(dist2==min(dist2),:) + y(dist2==min(dist2),:))/2; % midpoint
            dist2min(t,:) = [min(dist2) midCoord(1,:)]; % min of single voxels
        end

        % get mean coordinate of minimum distance
        dist2found = dist2min(dist2min(:,1)==min(dist2min(:,1)),:);
        if isvector(dist2found) % only one minimum of miminums
            dist2found = round(dist2found);
        else % might be more than one minimum of minimums
            dist2found = round(mean(dist2found));
        end

        % storing [nearby stone ID, distance, contact and centroid coordinates
        partCent = round(dataParticles(tab(j,1)).Centroid);
        contact(j,:) = [tab(j,1) dist2found partCent];

    end

    fprintf('...writing contact points for ParticleID = %u\n',k)
    dataParticles(k).Contacts = contact;

end

save(outFile,'dataParticles','-append')

poolobj = gcp('nocreate');
delete(poolobj);

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
