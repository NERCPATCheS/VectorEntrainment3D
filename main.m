%%  Main script for processing images of stones and their entrainment metrics
% This is the main script for running the VectorEntrainment3D suite of
% subroutines and functions used to calculate grain characteristics,
% their spatial orientation and proximity, and the entrainment of surface
% grains and their associated metrics from X-ray computed tomography
% (XCT) scanned images of water-worked sediment at the channel surface.
% Two preliminary scripts are necessary (STEP 1 and STEP 2 below) to build
% two MAT files that are later used by the main script (STEP 3).
%
% Subroutines below rely on a path to processed images, a path for output
% images, and their corresponding MAT filenames to store image data and
% entrainment metrics.  The user must create these folders and build a
% 'filepath' MAT file.  The input parent folder is the '/stacks/' folder
% while the output parent folder is the '/out/' folder.  Please create
% your own 'sample' subfolders (e.g. /R1/ and /R2/) for each XCT scanned
% sample used in your project within each of these parent folders.
%
% Processed entrainment metrics are derive from two spatially matched images.
% Each sample folder located in the '/stacks/' folder must contain two
% subfolders, '/ParticlesFull/' and '/Matrix/', each of which contains
% TIFF image stacks of 2D overhead view slices of the bed of grains. The
% first image must be an 8-bit binary TIFF stack of separated particles
% (separated using imaging software) located in the '/ParticlesFull/'
% subfolder.  The second must be a 16-bit grayscale TIFF stack of the
% corresponding fine-grained matrix that was separated from the particles
% during image segmentation (using appropriate imaging software). This
% image set must be located in '/Matrix/' subfolder.  It is important that
% the two stack image sets have the same spatial position, resolution and
% dimensions (please see stack images for R1 and R2 as examples).
%
% Be sure the image reference frame has a direction for water flow
% along the positive x-axis and gravity along the positive z-axis.
% Image must be a top down view such that layers are normal to the z-axis.
% During the STEP3 run, the code creates two more folders, 'Rotated' and
% 'Surface', in the '/out/' parent folder for output processed images.
%
% Other parameter requirements are listed below in the STEP 3 script.
% Please see details in the README.md file located on the PATCheS Project
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
%
% ----------------------------------------------------------------------
% STEP 3 processing order and description of entrainment model scripts:
%
%     (1) ImgStacks: Imports TIFF image stacks into the workspace
%
%     (2) ImgContacts: Searches contact points between stones
%
%     (3) ImgParticles: Calculates granular metrics for each stone
%
%     (4) ImgBedExtend: Extends bed of sample to reduce edge effects
%
%     (5) ImgSurfaces: Creates rotated images for exposure calculations*
%
%     (6) ImgExposure: Calculates exposure factor from ImgSurfaces output
%
%     (7) ImgEntrainment: Calculates entrainment metrics for each grain
%
% *NOTE: ImgSurfaces and ImgExposure must be run sequentially as return
% objects from the former are used as subroutine arguments for the latter.
%
% Please see details in the README.md file located on the PATCheS Project
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
% ----------------------------------------------------------------------
% LICENSING
% Copyright (C) 2018  PATCheS Project (http://www.nercpatches.org/)
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.  Please cite
% Voepel et al (2019) below for code use or modification and publication.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% REFERENCES
% Voepel, H., J. Leyland, R. Hodge, S. Ahmed, and D. Sear (2019),
% Development of a vector-based 3D grain entrainment model with
% application to X-ray computed tomography (XCT)scanned riverbed
% sediment, Earth Surface Processes and Landforms, doi: 10.1002/esp.4608
%
% Schmeeckle, M. W., J. M. Nelson, and R. L. Shreve (2007), Forces on
% stationary particles in near-bed turbulent flows, J. Geophys. Res.,
% 112, F02003, doi:10.1029/2006JF000536.

%% STEP 1: build filepaths and MAT filename (RUN BEFORE MAIN CODE)

clc
clear

filePathXlsx = char(strcat(pwd,'/filepaths.xlsx'));
[num,txt,~] = xlsread(filePathXlsx,'data');
[~,header,~] = xlsread(filePathXlsx,'header');

stacksPath = txt(:,1);
outPath	= txt(:,2);
outFile	= txt(:,3);
idxBottom = num(:,1);

clear filePathXlsx num txt
save filepath

%% STEP 2: build sieving database for sample GSD (RUN BEFORE MAIN CODE)

clc
clear

sizes = xlsread('SurfTableData.xlsx','Sizes');
[~,headers,~] = xlsread('SurfTableData.xlsx','Headers');
data = xlsread('SurfTableData.xlsx','SurfData');
[~,textDat,~] = xlsread('SurfTableData.xlsx','SurfInfo');

save SurfTableData

%% STEP 3: processing images and calculating entrainment (MAIN CODE)


clc
close all
clear

load filepath

fileNum = [1 2]; % samples [R1 R2]

% parameters for contact search algorithm
seDiam = 15; % diameter of spherical structural element
sampleRate = 0.1; % size of resample of object stone
maxDist2 = 30; % maximum squared voxel distance between stones 25
outGroup = 0; % = 0 writes nothing to file

% grain size percentiles required for subroutines (see STEP 2 above)
[~, ~, d50, ~, d84] = gsdData('SurfTableData', 999, false);

% parameters for entrainment model
resMM = 0.5999; % resolution (mm) of image
radTrim = round(d50/resMM); % (d50 as index) trim to overlap stones
Cd = 0.91; % coefficient of drag force (Schmeeckle et al, 2007)
Cl = 0.20; % coefficient of lift force
tol = 1e-5; % tolerance for zero moment sum
coheFlag = [false true]; % flag to use cohesive force (true applies Fc)
kill = 100; % max consecutive buried stone count before killing subroutine

% make viewing angles for 2D surface images
viewAngles = 0:15:90;

% setting overall timer
timerOverall = tic;

for k = fileNum

    timerVal = tic; % setting loop timer

    disp('************************************************************')
    fprintf('             Calculations for %s\n',outFile{k})
    disp('************************************************************')

    % STEP 1. Import binary images from stacks and create label matrices
    disp(['Importing binary images from image stacks for ' outFile{k} '...'])
    ImgStacks(stacksPath{k}, outPath{k}, outFile{k}, idxBottom(k), resMM)

    % STEP 2. Get contact particle IDs, distances, contact and grain centroids
    disp(['Mapping contact points between particles for ' outFile{k} '...'])
    ImgContacts(outFile{k}, seDiam, sampleRate, outGroup)

    % STEP 3. Get metrics of particle characteristics for each whole stone
    disp(['Getting metrics of particle characteristics for ' outFile{k} '...'])
    ImgParticles(outFile{k}, resMM)

    % STEP 4. Generate random wall within particle label matrix
    disp(['Adding extended bed to full particle label matrix for ' outFile{k} '...'])
    ImgBedExtend(outFile{k}, outPath{k}, radTrim);

    % STEP 5. Construct 2D surfaces of label particles at viewing angles
    disp(['Making 2D surface views of particle label matrix for ' outFile{k} '...'])
    [viewImages, surfaces] = ImgSurfaces(outFile{k}, outPath{k}, viewAngles);

    % STEP 6. Construct exposure factors EF for each stone [REQUIRES STEP 5 OUTPUTS]
    disp(['Running ImgExposure for ' outFile{k} '...'])
    ImgExposure(outFile{k}, viewAngles, viewImages, surfaces, resMM);

    % STEP 7. Running pivot angles separately to calibrate with tilt table
    disp(['Running ImgEntrainment for ' outFile{k} '...'])
    ImgEntrainment(outFile{k}, maxDist2, d50, d84, Cd, Cl, resMM, tol, kill, coheFlag(k))

    fprintf('\n\n')

    disp('************************************************************')
    fprintf('       %s completed on %s\n', outFile{k}, char(datetime))
    fprintf('           Database lapse time is %0.2f minutes\n',toc(timerVal)/60)
    disp('************************************************************')

    fprintf('\n\n\n\n')

end

disp('************************************************************')
fprintf('          Overall lapse time is %0.2f minutes\n',toc(timerOverall)/60)
disp('************************************************************')
