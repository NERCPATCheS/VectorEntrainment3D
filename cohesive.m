% Empirical cohesive force model

% Experiment was factorial design with response of vertical force (Y) pulls
% of a 25mm glass marble from a wet sand-clay mixuture where the matrix had
% sufficient time to drain to soil retention.  Experiment was designed with
% 100% clay as a control to evaluate two treatment effects: sand size and
% fraction of clay in the mixture (see Appendix C of Voepel et al, 2019). 

% We develop a power law model that calculates force per stone surface area 
% given sand size and fraction of clay content in the matrix.  We need to
% divide force/surface area of mixture by mean(force/surface area) of clay.
% First check whether force pulls from 100% clay are linearly dependent on 
% burial depth.  If so, regress three power law models, one for each burial
% depth using 3 corresponding mean(forceClay/area) divisors for each depth.  
% If not, regress single model using common mean(forceClay/area) divisor.

% Exponentiate linearised model to return to power law remembering to
% multipy exp(intercept) by mean(forceClay/area) to obtain model slope.

clc
clear
close all

load cohesive


% STEP 1: CHECK DEPTH DEPENDENCY OF FORCE PULLS FOR CLAY
% (NOT SIGNIFICANT, pVal = .526)
y = clay(:,2);
n = length(y);
X = [ones(n,1) clay(:,1)];
[~,~,~,~,statsClay] = regress(y,X); % stats = [R2, F, pval, sig2]

% Force pulls for clay not dependent on burial depth, using a single model.
% For a linearised power law log(Y) = log(a) + b*log(D) + c*log(F) + error,
% we have:
%           log(Y) = log(mixForce/area) - log(mean(clayForce/area))
%           log(D) = log(mean size for range of sand sizes)
%           log(F) = log(clay fraction used in mixture)
%            error ~ N(0,sig2)  CHECK THIS ASSUMPTION

% STEP 2: LINEAERISED POWER LAW MODEL REGRESSION
% need to adjust zero fraction values by adding 0.001
frac = mix(:,3);
frac(frac == 0) = 0.001;
mix(:,3) = frac;

% linear model
% (NOT SIGNIFICANT, pVal < .001, R2 = 0.54)
y = log(mix(:,7));
n = length(y);
X = [ones(n,1) log(mix(:,2:3))];
[b,bint,r,rint,stats] = regress(y,X); % stats = [R2, F, pval, sig2]

% checking residuals for normality (looks normal!)
% qqplot(r)

% STEP 3: PLOTTING PAIRED COMPARISON
% power law model intercept
a = exp(b(1))*mean(clay(:,4));

% get unique values for size and frac
size = unique(mix(:,2));
frac = unique(mix(:,3));

% generate mean Y and Yhat values to compare
muHat = zeros(12,4);
idx = 0;
for ii = 1:3
    for jj = 1:4
        idx = idx + 1;
        muHat(idx,1) = size(ii);
        muHat(idx,2) = frac(jj);
        muHat(idx,3) = mean(mix(mix(:,2)==size(ii)&mix(:,3)==frac(jj),6));
        muHat(idx,4) = a*size(ii)^b(2)*frac(jj)^b(3);
    end
end

% generate paired comparison plot
loglog([0.0001, 0.01],[0.0001, 0.01],'-k')
hold on
h = loglog(muHat(:,3),muHat(:,4),'Marker','o');
h.Color = hex2rgb('#d95f02');
h.Marker = 'o';
h.MarkerFaceColor = hex2rgb('#d95f02');
h.MarkerSize = 8;
h.LineStyle = 'none';
hold off
xlim([0.0001, 0.01])
ylim([0.0001, 0.01])
str{1} = sprintf('a = %0.6f', a);
str{2} = sprintf('b = %0.4f', b(2));
str{3} = sprintf('c = %0.4f', b(3));
text(0.002,0.0003,str,'FontSize',20,'Color','black','HorizontalAlignment','Left')
xlabel('observed force/area (N mm^{-2})')
ylabel('fitted force/area (N mm^{-2})')
set(gca,'FontSize',20)


% writing output files (comment for graph display)
% saveas(gcf,'figC1','epsc')



%% import data and save to MAT files (comment when finished)

% clc
% clear
% close all
% 
% 
% 
% % read in data and headers
% [clay, clayHeader, ~] = xlsread('clay.xlsx');
% [mix, mixHeader, ~] = xlsread('mix.xlsx');
% 
% % calculating force per area from depth for both clay and mix
% clay = [clay 2*pi*(25/2)*clay(:,1)]; % area for spherical cap
% mix = [mix 2*pi*(25/2)*mix(:,1)];
% clay = [clay clay(:,2)./clay(:,3)]; % force/area
% mix = [mix mix(:,4)./mix(:,5)];
% mix = [mix mix(:,6)/mean(clay(:,4))]; % (force/area)/mean(clayforce/area)
% 
% % updating headers
% clayHeader{3} = 'areaMM2';
% mixHeader{5} = 'areaMM2';
% clayHeader{4} = 'forcePerArea';
% mixHeader{6} = 'forcePerArea';
% mixHeader{7} = 'forcePerAreaPerMeanClay';
% 
% % save data for use above
% save('cohesive.mat')


%%

clear
clc
load('Run3-4BB-ERsin2','dataParticles')

n = length(dataParticles);
tau = zeros(n,2);
for k = 1:n
    
    if(~isempty(dataParticles(k).Entrainment))
        tau(k,1) = k;
        tau(k,2) = dataParticles(k).Entrainment.TauCrStar;
    end
    
end

tau = tau(tau(:,1)>0,:);
