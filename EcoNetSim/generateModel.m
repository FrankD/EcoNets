function model = generateModel(simData,comData),
% This function generates the food web, dispersal matrix and growth rates
% for the community simulation


if exist(simData.filename.model,'file') | exist([simData.filename.model '.mat'],'file'),
    disp(['Initialising model from file: ' simData.filename.model])
    L = load(simData.filename.model);

    % Check to see if file contains the correct data
    if ~any(strcmp(fieldnames(L),'model')),
        error(['File ' simData.filename.model ' does not contain model definition'])
    elseif abs(size(model.initialX)-[comData.nSpecies, comData.nSites])>1e-10,
        error(['Model input file (' simData.filename.Model ...
            'does not have expected spatial extent or number of species.'])
    else
        model = L.model;
    end
else
    %% Growth rates of the species
    [model.speciesGrowthRates,model.resources] = growthRates(comData,simData);
    %% Growth rates of the species
    [model.speciesIntraDD,model.intraDD] = intraDDGam(comData,simData);

    %% Dispersal matrix
    disp('Generating dispersal matrix')
    model.dispersalMatrix = sparse(DispersalMatrix(comData.m,comData.nSites));

    %% Food web topology
    % based upon the niche model (Wiliams & Martinez 2000, Nature)
    disp('Generating food web topology')
    [model.webTopology, model.trophicLevel] = ...
        nicheModel(comData.nSpecies,comData.connectance,false,simData.RandomSeed.niche);
    % Now generate interaction strengths for this web.

    [model.speciesInteractionStrengths] = interactionStrengths(model.webTopology,comData,simData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in initial abundance data if it exists
if exist(simData.filename.initialX,'file') | exist([simData.filename.initialX '.mat'],'file'),
    disp(['Initialising initial abundance from file: ' simData.filename.initialX])
    
    L = load(simData.filename.initialX);

    % Check to see if number of species and spatial area is that expected
    if ~any(strcmp(fieldnames(L),'X')),
        error(['File ' simData.filename.initialX ' does not contain model definition']);
    elseif abs(size(L.X)-[comData.nSpecies, comData.nSites])>1e-10,
        error(['Model input file (' simData.filename.initalX ...
            'does not have expected spatial extent or number of species.'])
    else
        model.initialX = L.X;
    end

else
    %% Initial Abundance
    model.initialX = initialAbundance(comData,simData);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersal matrix function
function M = DispersalMatrix(m,nSites),
% Calculate dispersal matrix
% M(i,j) is the probability that an individual moves from site i to site j
type = 'constrained';

if m<eps,
    % if m is very small assume no dispersal
    M = eye(nSites);
else
    k = -log((sqrt(1+4*m^2)-1)/2/m);
    if strncmp(type,'inf',3)
        %  Calculate decay const for exponetial dispersal
        k = -log((sqrt(1+4*m^2)-1)/2/m);
        sitePos = repmat([1:nSites],nSites,1);
        M = (1-exp(-k))/(1+exp(-k))*exp(-k*abs(sitePos-sitePos'));
        M(M<1e-10) = 0;
    elseif strncmp(type,'con',3),
        % Migration is only allowed to happen within environment
        pos = [1:nSites];
        for n=1:nSites,
            M(n,:) = exp(-k*abs(n-pos)) / sum(exp(-k*abs(n-pos)));
        end
        M(M<1e-10) = 0;
        % Renormalise
        M = M ./ repmat(sum(M,2),1,nSites);
    end
end
return


%% Growth rate function
function [r,resources] = growthRates(comData,simData),

if simData.RandomSeed.growthRate~=0,
    % Set random seed if required
    s = rand('state');
    sn = randn('state');
    rand('state',simData.RandomSeed.growthRate(1));
    randn('state',simData.RandomSeed.growthRate(2:3)');
end

% Resource coordinates
X = linspace(0,1,comData.nSites);

% Pick the number of resources that each species depends upon
resources.nPerSpecies = ceil(rand(comData.nSpecies,1)*comData.nResources);
% Finally calculate the growth rate of the different species at each point
% in space
for i=1:comData.nSpecies,
    % Now pick the optimal values for these resources for each species
    resources.optimalValues(i) = rand(1,1);
    r(i,:) = comData.rMin + (comData.rMax-comData.rMin) * ...
        exp(-(X-resources.optimalValues(i)).^2 / comData.rSigma *resources.nPerSpecies(i));
end

% % Method of creating species growth rates
% % Method = 1  random, independent of space
% % Method = 2  product of gaussians
% % Method = 3  sum of linear contributions
% method = 3;
% 
% if method==1,
%     % Growth rates for each species independent of spatial position
%     r = comData.r_mean * abs(1+sqrt(comData.r_cv)*repmat(randn(comData.nSpecies,1),1,comData.nSites));
%     % Rescale the growth rates
%     r = r / mean(nonzeros(r));
% elseif method==2,
%     % For each resource define the position at which optimum resource occurs
%     optimalX =  0.5*sqrt(comData.nResources)*randn(comData.nSpecies,comData.nResources);
%     % weighting of each resource to the growth rate
%     weights = 10*rand(comData.nSpecies,comData.nResources);
% 
% 
%     % This lets growth rate be the product of N Gaussian distributions
%     r = exp(- (sum(weights,2) * X.^2 ...
%         - 2*sum(weights.*optimalX,2) * X + ...
%         sum(weights.*optimalX.^2,2)*ones(size(X)))/comData.nResources);
%     r_mean = mean(r,2);
%     r = comData.r_mean * r ./ r_mean(:,ones(1,comData.nSites));
% elseif method==3
%     % Let the growth rate be a sum of linear contributions
%     % First assumes that growth rate must be positive
%     r_x0 = mean(max([0 comData.rMin]) + (comData.rMax - max([0 comData.rMin])) * rand(comData.nSpecies,comData.nResources),2);
%     % Second assumes it can have the full range
%     r_x1 = mean(comData.rMin + (comData.rMax - comData.rMin) * rand(comData.nSpecies,comData.nResources),2);
% 
%     % Place positive growth rate on left (X=1)
%     resources.r_x0 = r_x0;
%     resources.r_x1 = r_x1;
% 
%     ind = rand(comData.nSpecies,1)<0.5;
%     % Half the time place positive growth rate on right (X=max(X))
%     resources.r_x0(ind,:) = r_x1(ind,:);
%     resources.r_x1(ind,:) = r_x0(ind,:);
% 
%     r = resources.r_x1(:,ones(size(X))) + ...
%         ( resources.r_x0 - resources.r_x1 ) * (X - min(X))/(max(X)-min(X));
% end


if simData.RandomSeed.growthRate~=0,
    rand('state',s);
    randn('state',sn);
end

return

%% Intra-DD gam
function [gam, intraDD] = intraDDGam(comData,simData),
% A function to give the itraspecific density dependence

if simData.RandomSeed.intraDD~=0,
    % Set random seed if required
    s = rand('state');
    rand('state',simData.RandomSeed.intraDD);
end


% Resource coordinates
X = [1:comData.nSites];

% Let the gam be a sum of linear contributions
% gam must be positive
gam_x0 = mean(comData.gamMin + (comData.gamMax - comData.gamMin) * rand(comData.nSpecies,comData.nResources),2);
gam_x1 = mean(comData.gamMin + (comData.gamMax - comData.gamMin) * rand(comData.nSpecies,comData.nResources),2);

intraDD.gam_x0 = gam_x0;
intraDD.gam_x1 = gam_x1;

gam = intraDD.gam_x1(:,ones(size(X))) + ...
    ( intraDD.gam_x0 - intraDD.gam_x1 ) * (X - min(X))/(max(X)-min(X));

if simData.RandomSeed.intraDD~=0,
    rand('state',s);
end


return

%% Interaction strengths
function I = interactionStrengths(webTopology,comData,simData)

if simData.RandomSeed.interactionStrengths~=0,
    % Set random seed if required
    s = randn('state');
    randn('state',simData.RandomSeed.interactionStrengths');
end

% Make them log-normally distributed

% Define logicals if interactions are zero
noInteractionNeg = abs(comData.interaction.NegMean)<1e-10;
noInteractionPos = abs(comData.interaction.PosMean)<1e-10;
if nnz(webTopology) & (~noInteractionNeg | ~noInteractionPos),

    I = -webTopology + webTopology';
    I =  spdiags(spdiags(webTopology,0),0,I);

    % Calculate mean and variance of normal distribution giving rise to
    % lognormal
    if noInteractionPos,
        I(I>0) = 0;
    else
        PosS2 =log(comData.interaction.PosSigma2 + comData.interaction.PosMean^2) - ...
            log(comData.interaction.PosMean^2);
        PosM = log(comData.interaction.PosMean) - PosS2/2;
        I(I>0) = exp(PosM + sqrt(PosS2)*randn(nnz(I>0),1));
    end

    if noInteractionNeg,
        I(I<0) = 0;
    else
        NegS2 =log(comData.interaction.NegSigma2 + comData.interaction.NegMean^2) - ...
            log(comData.interaction.NegMean^2);
        NegM = log(comData.interaction.NegMean) - NegS2/2;
        I(I<0) = -exp(NegM + sqrt(NegS2)*randn(nnz(I<0),1));
    end
else
    I = sparse(zeros(comData.nSpecies));
end

if simData.RandomSeed.interactionStrengths~=0,
    randn('state',s);
end

return

%% Initial Abundance
function initialX = initialAbundance(comData,simData)

if simData.RandomSeed.initialAbundance~=0,
    % Set random seed if required
    s = rand('state');
    rand('state',simData.RandomSeed.initialAbundance);
end


initialX = comData.X_Init * ones(comData.nSpecies,comData.nSites);
% Keep those species who have gone extinct
ind = initialX < comData.Xmin;
initialX(ind) = 0;

% % And remove species so that the community size on average is
% % comData.nSpecies_mean
% for k=1:comData.nSites,
%     numSp = ceil(comData.nSpecies_mean*(1+randn(1)));
%     if numSp>comData.nSpecies,
%         numSp=comData.nSpecies;
%     elseif numSp<0,
%         numSp = 0;
%     end
%     rand_sp = randperm(comData.nSpecies);
%     initialX(rand_sp(1:(comData.nSpecies-numSp)),k)= 0;
% end

% Make X and c sparse matrices
initialX = sparse(initialX);
return
