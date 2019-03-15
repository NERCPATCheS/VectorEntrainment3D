
function ImgEntrainment(outFile, maxDist2, d50, d84, Cd, Cl, resMM, tol, kill, coheFlag)
%
% ImgEntrainment calculates entrainment metrics for sample surface grains.
%
% ImgEntrainment(outFile, maxDist2, d50, d84, Cd, Cl, resMM, tol, kill, coheFlag)
% is a vector-based, moment balance entrainment model for calculating critical
% shear stress and associated metrics for coarse surface grains.  Entrainment
% is calculated from the moment balance about a 3D axis of rotation (AOR).
%
% ImgEntrainment uses a binary search algorithm to calculate the moment balance.
% Pivot angle and the two plane of rotation (POR) angles for bearing and tilt
% are calculated once critical shear stress is calculated.
%
% ImgEntrainment requires the MATLAB Parallel Computing Toolbox and calls
% on the following subroutine arguments:
%
%   outFile = name of MAT file to store the updated 'dataParticles' array
%   maxDist2 = maximum squared voxel distance allowed to contact grains
%   d84 = the 84th percetile of grain size distribution of sample grains
%   Cd = coefficient of hydraulic drag force for natural stones
%   Cl = coefficient of hydraulic lift force for natural stones
%   resMM = voxel side length resolution for the scanned sample image (mm)
%   tol = tolerance of moment balance used to obtain critical shear stress
%
% Each stone in the sample is tested for 'potential entrainment' by testing
% for nonzero entrainment force.  If entrainment is possible, ImgEntrainment
% creates a structure array 'Entrainment', which is attached to the structure
% array 'dataParticles', and contains the following metrics:
%
%   PivotAngle = angle within the POR between gravity and the AOR (degrees)
%   PlaneBearing = departure angle of POR from the downstream vector (degrees)
%   PlaneTilt = departure angle of POR from the gravity vector (degrees)
%   ProjectionMM = distance from mean local bed elevation to top of grain (mm)
%   TauCr = critical shear stress for grain entrainment (Pa)
%   TauCrStar = dimensionless critical shear stress (-)
%   DragForce = hydrodynamic fluid drag force (N)
%   LiftForce = hydrodynamic fluid lift force (N)
%   CohesiveForce = tensile resistance force due to fine-grain matrix (N)
%   SubmergedForce = submerged weight force of grain (N)
%   LeftRightPartIDs = left and right contact pairs forming rotation axis
%   Tests = proportions of four viable contact pair tests plus overall test
%   Params = a Matlab structure array of entrainment parameter settings
%
% Entrainment parameter settings in the 'Params' structure are:
%
%   MaxDist2 = maximum squared distance between object and contact stones
%   D50MM = the 50th percentile of surface grain size (mm)
%   D84MM = the 84th percentile of surface grain size (mm)
%   Cd = drag force coefficient for drag force calculations
%   Cl = lift force coefficient for lift force calculations
%   ResMM = image voxel resolution (mm)
%   Tol = moment balance convergence tolerance
%   Kill = consecutive non-entrained grain count before killing algorithm
%
% Please see details in the README.md file located on the PATCheS Project
% GitHub page (https://github.com/NERCPATCheS/VectorEntrainment3D).
%
% AUTHOR: Hal Voepel
% DATE: 15 October 2018
%
% See also ImgStacks, ImgContacts, ImgParticles, ImgBedExtend, ImgSurfaces,
% and ImgExposure.

% REFERENCES
% Voepel, H., J. Leyland, R. Hodge, S. Ahmed, and D. Sear (2019),
% Development of a vector-based 3D grain entrainment model with
% application to X-ray computed tomography (XCT)scanned riverbed
% sediment, Earth Surface Processes and Landforms, doi: 10.1002/esp.4608
%
% Copyright (C) 2018  PATCheS Project (http://www.nercpatches.org/)


%---------CHECKING REQUIREMENTS BEFORE RUN------------

% Check user has required toolbox(s) installed installed
tbCheck = license('test', 'Distrib_Computing_Toolbox');
if ~tbCheck
	% User does not have the toolbox installed.
	error('Requires Parellel Computing Toolbox.')
end

% loading database
load(outFile);

% getting number of particles
partCount = ccParticlesFull.NumObjects;

% initiallizing loop break counter
killCounter = 0;


%-----------REMOVING PREVIOUS STRUCURE FROM DATAPARTICLES----------
try
    dataParticles = rmfield(dataParticles,'Entrainment');
    disp('Removing previous Entrainment field from dataParticles')
catch ME
end


%----------------DEFINE PHYSICAL CONSTANTS--------------------
g = 9.81; % gravitational acceleration (m/s2)
rho = 1000; % density of water (kg/m3)
rhoS = 2650; % density of sediment (kg/m3)
cohesiveArea = GetCohesiveArea(999); % cohesive force per area (N/mm2)

%----------PREALLOCATING STRUCTURES FOR DATAPARTICLES-----------------
Entrainment = struct('PivotAngle',-1,'PlaneBearing',-1,'PlaneTilt',-1,'ProjectionMM',-1,...
    'TauCr',-1,'TauCrStar',-1,'DragForce',-1,'LiftForce',-1,'CohesiveForce',-1,...
    'SubmergedWeight',-1,'LeftRightPartIDs',-1,'Tests',-1,'Params',-1);

Params = struct('MaxDist2',-1,'D50MM',-1,'D84MM',-1,'Cd',-1,'Cl',-1,'ResMM',-1,'Tol',-1, 'Kill',-1);


%----------PASSING DATA TO PARALLEL POOL WORKERS-----------
numCores = feature('numcores');
poolobj = parpool(numCores);
data = load(outFile); % passing array to workers
labelParticlesBed = data.labelParticlesBed;
dataParticles = data.dataParticles;


for k = 1:partCount % looping over particles

    % starting lapse timer
    tic

    %-----------DROP STONES THAT ARE BOTTOM CROPPED------------
    if (dataParticles(k).CroppedFlag)

        disp('=========================================================')
        fprintf('Stone bottom cropped for ParticleID = %u\n', k)
        killCounter = killCounter + 1; % increment kill switch
        continue

    end

    %-----------CHECKING STONE EXPOSURE TO FORCES--------------
    [Fd, Fl, projMM, idxZ] = EntrainForces(dataParticles(k), ...
        labelParticlesBed, d84, 1000, Cd, Cl, 0, 0, resMM);

    % If no force, then continue to next particle
    if (sum(sum([Fd Fl])) == 0 || isnan(idxZ))

        disp('=========================================================')
        fprintf('No entrainment forces for ParticleID = %u\n', k)

        % kill loop when more than 20 consecutive particles are buried
        killCounter = killCounter + 1;
        if killCounter > kill
            fprintf('More than %u consecutive buried stones...\n', killCounter)
            fprintf('killing process for remaining stones in %s\n', outFile)
            break
        end
        continue

    end


    %-----------GETTING VIABLE CONTACT LIST--------------------
    [pairs, tests] = GetContactPairs(dataParticles(k),maxDist2);

    % checking for nonexistant or single contact point
    if (isscalar(pairs) || isempty(pairs))

        disp('=========================================================')
        fprintf('No contact points for ParticleID = %u\n', k)
        killCounter = killCounter + 1; % increment kill switch
        continue

    else % particle found with entrainment

        % getting pairs count reset kill loop counter
        disp('=========================================================')
        fprintf('Entrainment structure created for ParticleID = %u\n', k)
        n = size(pairs,1); % there are at least one pair
        killCounter = 0; % resetting kill switch

    end


    %-----------CALCULATING RESISTIVE FORCE---------------------
    % getting metrics from database
    vol = dataParticles(k).VolumeMM3/1000^3; % particle volume (m3)
    D = dataParticles(k).MedianAxisMM/1000; % median axis diameter (m)
    cent = dataParticles(k).Centroid; % particle centroid
    mat = dataParticles(k).MatrixContactAreaRatio; % matrix cont area ratio
    bwSurf = labelParticlesSurf == k; % surface element count for kth particle
    surfArea = sum(bwSurf(:))*resMM^2; % surface area of stone (mm2)

    % calculate forces and sum results as resistive force
    if coheFlag
        Fc = [0 0 cohesiveArea*surfArea*mat]'; % cohesive force (N), EQ (14a)
    else
        Fc = [0 0 0]'; % do not apply cohesive force
    end
    Fw = [0 0 (rhoS - rho)*g*vol]'; % submerged weight of water (N), EQ (13)
    Fr = Fc + Fw; % resistive force (N)


    %-------------CALCULATING ENTRAINMENT MOMENT ARM-----------------
    Re = zeros(n,3); % preallocating for entrainment arm
    cent = [cent(1:2) idxZ]; % centroid of applied Fe point
    cont = dataParticles(k).Contacts; % getting table of contacts
    for j = 1:n
        Re(j,:) = cont(cont(:,1) == pairs(j,1),3:5); % left contact centroid
    end
    Re = repelem(cent,n,1) - Re; % moment arms for entrainment forces



    %----------GETTING CRITICAL SHEAR FOR EACH CONTACT PAIR----------
    fprintf('Finding critical shear from %u viable contact pairs\n',n)
    shear = zeros(n,1);

    parfor j = 1:n % looping through contact pairs


        %-------------RETRIEVING INFORMATION FOR CONTACT PAIRS----------
        % get rotation axis unit vector
        lambda = pairs(j,12:14)';

        % resistance moment arm from contact proj onto POR
        Rr = -pairs(j,18:20)';

        % get bearing and tilt angles for plane of rotation
        beta = pairs(j,22); % bearing angle
        gamma = pairs(j,23); % tilt angle


        %--------GETTING INITIAL ENTRAINMENT CALCULATIONS-------------
        % setting initial range for shear from no flow conditions
        dTau = 200; % N/m2
        tau = [0 dTau];

        % getting lower bound forces (N) for no flow (zero for both)
        [Fd1, Fl1, ~, ~] = EntrainForces(dataParticles(k), ...
            labelParticlesBed, d84, tau(1), Cd, Cl, gamma, beta, resMM);

        % getting upper bound forces (N) for dTau
        [Fd2, Fl2, ~, ~] = EntrainForces(dataParticles(k), ...
            labelParticlesBed, d84, tau(2), Cd, Cl, gamma, beta, resMM);

        % summing entrainment forces at boundaries
        Fe1 = Fd1 + Fl1; % lower bound of force
        Fe2 = Fd2 + Fl2; % upper bound of force

        % calculating initial moment at tau = [0 dTau], EQ (7)
        m01 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe1));
        m02 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe2));


        %---------RESETTING SHEAR INTERVAL TO ZERO MOMENT INTERVAL-------

        % setting iteration indices
        idx = 1;

        % adjusting tau value towards m0 = 0 results
        adj = floor(m01/(m01 - m02)); % adjusting multiplier
        tau = [adj - 1 adj]*dTau; % tau value near m0 = 0

        % skipping impossible convergence to m0 = 0
        if adj < 0 % condition for moving away from m0 = 0
            shear(j) = inf;
            continue
        end

        % getting lower bound forces (N) for no flow (zero for both)
        [Fd1, Fl1, ~, ~] = EntrainForces(dataParticles(k), ...
            labelParticlesBed, d84, tau(1), Cd, Cl, gamma, beta, resMM);

        % getting upper bound forces (N) for dTau
        [Fd2, Fl2, ~, ~] = EntrainForces(dataParticles(k), ...
            labelParticlesBed, d84, tau(2), Cd, Cl, gamma, beta, resMM);

        % summing entrainment forces at boundaries
        Fe1 = Fd1 + Fl1; % lower bound of force
        Fe2 = Fd2 + Fl2; % upper bound of force

        % calculating initial moment at tau = [0 dTau], EQ (7)
        m01 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe1));
        m02 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe2));
        m0 = [m01 m02]; % initial moment interval to check


        %-------------BINARY SEARCH FOR MOMENT BALANCE------------
        while (max(abs(m0)) > tol && idx < 10000)

            % adjusting tau interval until moment interval contains zero
            while (prod(m0) > 0)

                % incrementing values
                idx = idx + 1;
                tau = tau + dTau;

                % getting new bounds of entrainment forces
                % new lower bound
                [Fd1, Fl1, ~, ~] = EntrainForces(dataParticles(k), ...
                    labelParticlesBed, d84, tau(1), Cd, Cl, gamma, beta, resMM);

                % new upper bound
                [Fd2, Fl2, ~, ~] = EntrainForces(dataParticles(k), ...
                    labelParticlesBed, d84, tau(2), Cd, Cl, gamma, beta, resMM);

                % recalculating moments
                Fe1 = Fd1 + Fl1; % lower bound of force
                Fe2 = Fd2 + Fl2; % upper bound of force
                m01 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe1));
                m02 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe2));
                m0 = [m01 m02]; % initial moment interval to check, EQ (7)

            end % ending inner while loop

            % selecting which half of m0 interval contains zero
            if prod([m0(1) mean(m0)]) < 0

                % first half of increment contains zero
                tau = [tau(1) mean(tau)];

            else

                % second half of increment contains zero
                tau = [mean(tau) tau(2)];

            end

            % incrementing index
            idx = idx + 1;

            % getting new bounds of entrainment forces
            % new lower bound
            [Fd1, Fl1, ~, ~] = EntrainForces(dataParticles(k), ...
                labelParticlesBed, d84, tau(1), Cd, Cl, gamma, beta, resMM);
            % new upper bound
            [Fd2, Fl2, ~, ~] = EntrainForces(dataParticles(k), ...
                labelParticlesBed, d84, tau(2), Cd, Cl, gamma, beta, resMM);
            Fe1 = Fd1 + Fl1; % lower bound of force
            Fe2 = Fd2 + Fl2; % upper bound of force
            m01 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe1));
            m02 = lambda'*(cross(Rr,Fr) + cross(Re(j,:)',Fe2));
            m0 = [m01 m02]; % subsequent moment interval to check, EQ (7)

        end % ending outer while loop

        % storing shear value for jth contact pair
        shear(j) = mean(tau);

    end % ending for j loop over contact pairs

    %---------GETTING MINIMUM CRITICAL SHEAR FROM PAIRS------------
    pairs = [pairs shear];
    idxShear = find(shear == min(shear),1);
    finalPair = pairs(idxShear,[1:2 21:24]);

    %--------STORING FINAL CALCULATIONS TO DATABASE---------------
    contIDs = finalPair(1:2);
    alpha = finalPair(3);
    beta = finalPair(4);
    gamma = finalPair(5);
    tauCr = finalPair(6);
    tauCrStar = tauCr/((rhoS - rho)*g*D);

    fprintf('\nFor left-right contact pair, Cl = %u and Cr = %u,\n',contIDs(1),contIDs(2))
    fprintf('critical shear, tauCr* = %0.4f, for ParticleID = %u\n\n', tauCrStar, k)

    % getting final entrainment force values
    [Fd, Fl, ~, ~] = EntrainForces(dataParticles(k), ...
        labelParticlesBed, d84, tauCr, Cd, Cl, gamma, beta, resMM);

    % storing information into structures
    dataParticles(k).Entrainment.PivotAngle = alpha;
    dataParticles(k).Entrainment.PlaneBearing = beta;
    dataParticles(k).Entrainment.PlaneTilt = gamma;
    dataParticles(k).Entrainment.ProjectionMM = projMM;
    dataParticles(k).Entrainment.TauCr = tauCr;
    dataParticles(k).Entrainment.TauCrStar = tauCrStar;
    dataParticles(k).Entrainment.DragForce = Fd;
    dataParticles(k).Entrainment.LiftForce = Fl;
    dataParticles(k).Entrainment.CohesiveForce = Fc;
    dataParticles(k).Entrainment.SubmergedWieght = Fw;
    dataParticles(k).Entrainment.LeftRightPartIDs = contIDs;
    dataParticles(k).Entrainment.Tests = tests;
    dataParticles(k).Entrainment.Params.MaxDist2 = maxDist2;
    dataParticles(k).Entrainment.Params.D50MM = d50;
    dataParticles(k).Entrainment.Params.D84MM = d84;
    dataParticles(k).Entrainment.Params.Cd = Cd;
    dataParticles(k).Entrainment.Params.Cl = Cl;
    dataParticles(k).Entrainment.Params.ResMM = resMM;
    dataParticles(k).Entrainment.Params.Tol = tol;
    dataParticles(k).Entrainment.Params.Kill = kill;

    % displaying elapsed time for kth particle
    toc


end % ending loop over particles


% shutting down workers
disp('=========================================================')
delete(poolobj);


% saving dataParticles
save(outFile,'dataParticles','-append')
% savedFile = strcat(strtok(outFile,'.'),'linear.mat');
% save(savedFile,'dataParticles')

end % end function

%*********************LOCAL FUNCTION BLOCK**********************

% GetCohesiveArea

function [cohesiveArea] = GetCohesiveArea(fileNum)

%---------------CALCULATING COHESIVE FORCE PER AREA------------
% use overall GSD from sieving data to make weighted mean force/area

load('SurfTableData') % load sieving data
% clay (size < 63 um)       data(29)
% size (0.1 0.3] = 0.2      data(23:27)
% size (0.5 1.0] = 0.75     data(20:22)
% size (1.0 2.0] = 1.5      data(18:19)
sizes = [0.2 0.75 1.5]';
ss = sum(data(data(:,30) == fileNum, :))';
prop = [sum(ss(23:27)) sum(ss(20:22)) sum(ss(18:19))]';
frac = [prop repelem(ss(29),3,1)];
frac = frac./[sum(frac,2) sum(frac,2)];
frac = frac(:,2);
prop = prop/(sum(prop));

% tensile force per area (N/mm2) for binary size and fraction of clay, EQ (14b)
forceArea = @ (size,frac) 0.002579.*size.^-0.5332.*frac.^0.2328;

% weighted mean cohesive force per area (N/mm2), EQ (14c)
cohesiveArea = forceArea(sizes,frac)'*prop/sum(prop);

end



%EntrainForces

function [Fd, Fl, projMM, idxZ] = EntrainForces(dataPartK, lblBed, d84, tau, Cd, Cl, gamma, beta, resMM)

% converting d84 from mm to index
d84Idx = round(d84/resMM);

%------------------MAKING LOGARITHMIC ELEVATION PROFILE----------------

% gathering metrics for kth particle
bb = dataPartK.BoundingBox; % 3D box around particle
img = dataPartK.FrontalImage; % 2D downstream view of single stone
proTrue = round(dataPartK.MaxElevMM/resMM); % stone max elev as index
ohArea = dataPartK.OverheadAreaMM2*1e-6; % overhead area (m2) for lift force
er = dataPartK.ExposureRatio(1); % viewable fraction of ohArea
EF = dataPartK.ExposureFactor; % exposure factor (hiding effects), EQ (10)

% determining size of local bed
y = size(lblBed,2);
rngY = round(max(1,bb(1) - d84Idx)):round(min(y,bb(1) + bb(4) + d84Idx));
rngX = round(bb(2)):round(bb(2) + bb(5));

% clipping local bed z-vectors
localBed = lblBed > 0;
localBed = localBed(rngX,rngY,:);

% preallocating space for local bed elevations at (X,Y) coords
[x, y, ~] = size(localBed);
elevZ = zeros(x*y,1);

idx = 1;

for i = 1:x
    for j = 1:y
        try
            elevZ(idx) = find(localBed(i,j,:) > 0, 1);
            idx = idx + 1;
        catch ME
            idx = idx + 1;
        end
    end
end

% calculating mean bed
muBed = round(mean(elevZ)); % mean bed elevation as index
topElev = size(img,2); % highest elevation in image

% clipping frontal image to mean bed elevation
imgProt = img(:,(topElev-muBed):topElev); % projection above mean bed

% contruct velocity profile over projection image
z0 = d84Idx/10; % setting length scale to 10% of D84
f = log(((0:topElev)' + z0)/z0); % logarithmic profile
f = f(1:length((topElev-muBed):topElev));
% g = (0:(length(f)-1))*max(f)/(length(f)-1); % linear profile (for testing)


%---------------MAKING VELOCITY PROFILE AS FUNCTION OF SHEAR-----------

rho = 1000; % kg/m3
kappa = 0.407;
u2 = tau/rho/kappa^2*f.^2; % squared velocity profile, EQ (8)


% --------------------INTEGRATION OVER VELOCITY PROFILE--------------------
% THREE CASES: For p = projection height, D = stone bottom, z = 0 at mean bed
% 1. p > p - D > z = 0: stone above mean bed elevation (full force)
% 2. p > z = 0 > p - D: stone intersects mean bed elevation (partial force)
% 3. z = 0 > p > p - D: stone below mean bed elevation (no force)

% preallocating Riemann sum components
y = size(imgProt,2); % integration length to top of image
u2dA = zeros(y,1);

% initialize indices used over Riemann summation
idx = 1; % velocity profile index at z = 0
pmD = 1; % stone bottom index (to be adjusted)
p = 1; % protrusion index (to be adjusted)

if sum(img(:)) == sum(imgProt(:)) % 1. stone above mean bed elevation

    % finding idx distance from z = 0 to bottom of stone D
    while sum(imgProt(:,idx)) == 0

        idx = idx + 1;
        pmD = idx;

    end

    % reached bottom of stone at distance D, now looping over stone
    while sum(imgProt(:,idx)) > 0

        u2dA(idx) = u2(idx)*sum(imgProt(:,idx))*(resMM/1000)^2;
        idx = idx + 1;

    end

    p = idx; % top of protrusion

elseif sum(imgProt(:)) > 0 % 2. stone intersects mean bed elevation

    % storing rectangular Riemann sum components
    while sum(imgProt(:,idx)) > 0

        u2dA(idx) = u2(idx)*sum(imgProt(:,idx))*(resMM/1000)^2;
        idx = idx + 1;

    end

    p = idx; % top of protrusion

end % 3. stone below mean bed elevation if neither case above


% trimming vectors
h = 1:idx;
zElev = h; % protrusion height index
u2dA = u2dA(h);
u2 = u2(h);
f = f(h);

% calculating forces in vector form
Fu = sum(imgProt(:))/sum(img(:)); % stone areal fraction where u > 0
EFadj = 1 - max([0, (Fu - EF)/Fu]); % adjust EF, EQ (9b)
Fd = [Cd/2*rho*sum(u2dA)*EFadj 0 0]'; % drag force, EQ (9a)
Fl = [0 0 -Cl/2*rho*ohArea*er*(u2(p) - u2(pmD))]'; % lift force, EQ (11)

appFd = round(zElev*u2dA/sum(u2dA)); % profile distance to apply force, EQ (12)
idxZ = proTrue + (idx - appFd); % index distance from top to apply force
projMM = sum(u2dA > 0)*resMM; % projection height of stone (mm)



end

%GetContactPairs


function [pairs, tests] = GetContactPairs(dataPartK,maxDist2)


% checking existance of contacts
if isempty(dataPartK.Contacts)
    pairs = -1;
    tests = -1;
    return
end


% removing distance grains and getting particleIDs
cont = dataPartK.Contacts;
cont = cont(cont(:,2) <= maxDist2,:);

% check for less than 2 contacts
if (isempty(cont) || isvector(cont))
    pairs = -1;
    tests = -1;
    return
end
contID = cont(:,1);

% make contacts matrix and get size
cont = cont(:,3:5);
n = size(cont,1);

% particle and contact centoids
part = dataPartK.Centroid;
partMat = repmat(part,n,1);

% get contact angles from x-axis and z-axis
cp = cont - partMat; % center of mass to contact vector
phiX = atan2(cp(:,2),cp(:,1))*180/pi;
thetaZ = acos(cp(:,3)./sqrt(diag(cp*cp')))*180/pi;

% populate matrix with needed centroid information
cent = zeros(n,6);

cent(:,1) = contID; % contact stone ID
cent(:,2) = phiX; % angle from x-axis
cent(:,3) = thetaZ; % angle from z-axis
cent(:,4:6) = cont; % contact centroid

% allocating space for all possible contact pairs
n2 = nchoosek(n,2);
pairs = zeros(n2,26);

% allocating space for obstruction contact testing
% conQ = cell(n2,1); % MIGHT NOT NEED

% index for viableContacts
idx = 0;

% left-right orientation by sorting on PhiX
cent = sortrows(cent,2);

% build list of all left-right contact combinations
for i = 1:(n-1) % indexing left contact
    for j = (i+1):n % indexing right contact

        % increment tabling index
        idx = idx + 1;

        % getting left-right contact IDs
        pairs(idx,1) = cent(i,1); % left contact ID
        pairs(idx,2) = cent(j,1); % right contact ID

        % left-right contact vectors and their normal unit vector
        conL = cent(i,4:6) - part; % centering left contact on particle
        conR = cent(j,4:6) - part; % centering right contact on particle
        pairs(idx,3:5) = conL; % left contact vector
        pairs(idx,6:8) = conR; % right contact vector
        nCont = cross(conL,conR); % normal to contact plane
        pairs(idx,9:11) = nCont/norm(nCont); % normal unit vector

        % rotational axis and corresponding unit vector (lambda)
        rot = conL - conR; % rotation axis
        pairs(idx,12:14) = rot/norm(rot); % rotation axis unit vector

        % horizontal rotational axis and corresponding unit vector (lambdaXY)
        rotXY = [rot(1:2) 0]; % horizontal axis
        pairs(idx,15:17) = rotXY/norm(rotXY); % horizontal unit vector

        % orthogonal projection of contact vector onto plane of rotation
        projPOR = eye(3) - rot'*rot/(rot*rot');
        pairs(idx,18:20) = conL*projPOR;

        % pivot angle calculations (alpha), EQ (6a)
        projGrav = [0 0 9.81]*projPOR; % project gravity onto POR
        projCl = pairs(idx,18:20); % projection of Cl onto POR
        prodNorms = norm(projGrav)*norm(projCl); % product of norms
        pairs(idx,21) = real(acos(projGrav*projCl'/prodNorms)*180/pi);

        % bearing angle calculations (beta), EQ (6b)
        projPORXY = eye(3) - rotXY'*rotXY/(rotXY*rotXY');
        projX = [1 0 0]*projPORXY;
        pairs(idx,22) = real(acos([1 0 0]*projX'/norm(projX))*180/pi);

        % tilt angle calculations (gamma), EQ (6c)
        lambda =  pairs(idx,12:14);
        lambdaXY =  pairs(idx,15:17);
        pairs(idx,23) = real(acos(lambda*lambdaXY')*180/pi);

        % storing vector for TEST4: obstructing contact test
        pairs(idx,24:26) = -projX/norm(projX);

    end
end


%---------------TESTING PAIRS FOR VIABILITY------------------------

% TEST1: plane of rotation between contacts (left-right dot lambda)
leftVecDotLambda = dot(pairs(:,3:5),pairs(:,12:14),2);
rightVecDotLambda = dot(pairs(:,6:8),pairs(:,12:14),2);
test1 = leftVecDotLambda.*rightVecDotLambda < 0; % passes test

% TEST2: gravity topple over tilted POR (left-right dot lambdaXY)
leftVecDotLambdaXY = dot(pairs(:,3:5),pairs(:,15:17),2);
rightVecDotLambdaXY = dot(pairs(:,6:8),pairs(:,15:17),2);
test2 = leftVecDotLambdaXY.*rightVecDotLambdaXY < 0; % passes test

% TEST3: forward rotation | angle | < 90 degrees (proj'd contact dot [1 0 0])
xUnit = repmat([1 0 0],n2,1);
projCl = pairs(:,18:20);
test3 = dot(projCl,xUnit,2) > 0; % passes test

% TEST4: no obstructing contact in forward tilting area
test4 = false(n2,1); % preallocating logical vector

% checking each contact pair against all other contact points
for k = 1:n2

    % initiallize obstructing contact flag
    obsFlag = false;

    % getting table of all other contacts points aside from kth pair
    obstr = cent(cent(:,1) ~= pairs(k,1) | cent(:,1) ~= pairs(k,2),:);

    % replacing contact CM with vector centered on particle CM
    n3 = size(obstr,1);
    partMat = repmat(part,n3,1);
    obstr(:,4:6) = obstr(:,4:6) - partMat;

    % getting test vectors for kth contact pair
    nCont = pairs(k,9:11);
    nProj = pairs(k,24:26);

    % checking each contact point with kth pair
    for t = 1:n3
        obsVec = obstr(t,4:6);
        if obsVec*nCont' < 0 && obsVec*nProj' < 0
            obsFlag = true; % contact is obstructing pivot
            break
        end
    end

    if obsFlag
        test4(k) = false; % fails test
    else
        test4(k) = true; % passes test
    end

end


%----------------FINISHING LIST OF VIABLE PAIRS----------------------

% intersecting four tests
testAll = test1 & test2 & test3 & test4;

% storing testing results
tests = ([sum(test1) sum(test2) sum(test3) sum(test4) sum(testAll)]/n2)';

% generate viable list of contact pairs
% pairs = [pairs(testAll,:) pivotAngle(testAll) rotAngle(testAll) tiltAngle(testAll)];
% pairs = sortrows(pairs, [1 2]);
pairs = pairs(testAll,1:23);


end
