function model = generateModel(simData,comData),
% This function generates the food web, dispersal matrix and growth rates
% for the community simulation
%
% This model runs on a 2d spatial arena


if exist(simData.filename.model,'file') | exist([simData.filename.model '.mat'],'file'),
    disp(['Initialising model from file: ' simData.filename.model])
    L = load(simData.filename.model);

    % Check to see if file contains the correct data
    if ~any(strcmp(fieldnames(L),'model')),
        error(['File ' simData.filename.model ' does not contain model definition'])
    elseif abs(size(model.initialA)-[comData.nSpecies, comData.nSites])>1e-10,
        error(['Model input file (' simData.filename.Model ...
            'does not have expected spatial extent or number of species.'])
    else
        model = L.model;
    end
else
    disp('Generating dispersal matrix')
    if exist(simData.filename.DispersalMatrix,'file') | exist([simData.filename.DispersalMatrix '.mat'],'file'),
        disp(['Initialising dispersal matrix from file: ' simData.filename.DispersalMatrix])
        load(simData.filename.DispersalMatrix)
        model.dispersalMatrix = dispersalMatrix;
        model.xCoord = xCoord;
        model.yCoord = yCoord;
    else
        %% Dispersal matrix
        [model.dispersalMatrix, model.xCoord, model.yCoord] = DispersalMatrix(comData.m,comData.nSites, simData); 
    end
    %% Growth rates of the species
    %% (the final argument is the resource gradient)
    [model.speciesGrowthRates,model.resources] = growthRates(comData,simData,model.xCoord);
    %% Growth rates of the species
    [model.speciesIntraDD,model.intraDD] = intraDDGam(comData,simData,model.xCoord);


    %% Food web topology
    % based upon the niche model (Wiliams & Martinez 2000, Nature)
    disp('Generating food web topology')
    
    if(comData.connectance == -1)
      model.webTopology = zeros(4, 4);
%      model.webTopology(:, end) = ones(3, 1);
%      model.webTopology(end, end) = 0;
      model.webTopology(3, 1) = 1;
      model.webTopology(2, 3) = 1;
      model.webTopology(2, 1) = 1;
      model.trophicLevel = zeros(4, 1);
    else
      [model.webTopology, model.trophicLevel] = ...
       nicheModel(comData.nSpecies,comData.connectance,false,simData.RandomSeed.niche); % JY addition
%        nicheModel(comData.nSpecies,comData.connectance,false,simData.RandomSeed.niche,
%        simData); % JY edit: old code
    end


    % Now generate interaction strengths for this web.
    [model.speciesInteractionStrengths] = interactionStrengths(model.webTopology,comData,simData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in initial abundance data if it exists
if exist(simData.filename.initialA,'file') | exist([simData.filename.initialA '.mat'],'file'),
    disp(['Initialising initial abundance from file: ' simData.filename.initialA])

    L = load(simData.filename.initialA);

    % Check to see if number of species and spatial area is that expected
    if ~any(strcmp(fieldnames(L),'A')),
        error(['File ' simData.filename.initialA ' does not contain model definition']);
    elseif abs(size(L.A)-[comData.nSpecies, prod(comData.nSites)])>1e-10,
        error(['Model input file (' simData.filename.initalA ...
            'does not have expected spatial extent or number of species.'])
    else
        model.initialA = L.A;
    end

else
    disp('Calculating Initial Abundance')
    %% Initial Abundance
    model.initialA = initialAbundance(comData,simData);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dispersal matrix function
function [M,X,Y] = DispersalMatrix(m,nSites, simData),
% Calculate dispersal matrix
% M(i,j) is the probability that an individual moves from site i to site j
type = 'inf';

X = reshape(repmat([1:nSites(1)],nSites(2),1),prod(nSites),1);
Y = reshape(repmat([1:nSites(2)]',1,nSites(1)),prod(nSites),1);

if m<eps,
    % if m is very small assume no dispersal
    M = speye(prod(nSites));
else
    %  Calculate decay const for exponetial dispersal over a 2d arena
    % m is the mean dispersal distance
    nPoints = 1+2*max(nSites);
    x_tmp = repmat([1:nPoints],nPoints,1);
    d = ((reshape(x_tmp,nPoints^2,1)-(nPoints+1)/2).^2 + (reshape(x_tmp',nPoints^2,1)-(nPoints+1)/2).^2).^0.5;
    x0 = 2/m;
    k = fzero(@(x)m-sum(d.*exp(-x*d))/sum(exp(-x*d)),x0);
    %    k = 2/m; % The continuous result

    % Calculate normalisation constant as well
        N = sum(exp(-k*d));

    % Calculate squared distance between every site
    dist2 = (X(:,ones(size(X))) - X(:,ones(size(X)))').^2 + (Y(:,ones(size(Y))) - Y(:,ones(size(Y)))').^2;
    % Find separation distance that gives a prob of dispersal of 10*eps
    rMax = 2*log(k) - log(10^-10);
    % Convert distance matrix to a sparse matrix
    dist2(dist2>rMax^2) = 0;
    dist = sparse(dist2);
    dist = spfun(@sqrt,dist);

    if strncmp(type,'inf',3)
        M = spfun(@exp,-k*dist);
        dM = (M + speye(prod(nSites))) / N;
        % Again remove any entries that are too small
        %M(M<1e-10) = 0;
    elseif strncmp(type,'con',3),
        % Migration is only allowed to happen within environment
        for n=1:prod(nSites),
            if mod(n, round(prod(nSites)/100)) == 0 && simData.debug
                fprintf('Migration done %g\n', n/prod(nSites))
            end
            dM(n,:) = spfun(@exp,-k*dist(n,:));
            dM(n,n) = 1; % Make sure no dispersal can occur
            dM(n,:) = dM(n,:) / sum(dM(n,:)); % Renormalise to keep all individuals in the arena
        end
    end

    % Again remove any entries that are too small
      M = dM;
    %    M(M<1e-10) = 0;

    [ind1,ind2]  = find(M>1e-10);
    M = sparse(ind1,ind2,M(sub2ind(size(M),ind1,ind2)),prod(nSites),prod(nSites),length(ind1));
    % Renormalise after removing these entries
    M = M ./ repmat(sum(M,2),1,length(M));
end



return


%% Growth rate function
function [r,resources] = growthRates(comData,simData,R),

% R is the resource coordinate
% Make sure it is a row vector
R = reshape(R,1,length(R));

if any(simData.RandomSeed.growthRate~=0),
    % Set random seed if required
    s = rand('state');
    sn = randn('state');
    rand('state',simData.RandomSeed.growthRate(1));
    randn('state',simData.RandomSeed.growthRate(2:3)');
end


% Method of creating species growth rates
% Method = 1  random, independent of space
% Method = 2  product of gaussians
% Method = 3  sum of linear contributions
% Method 4 use 1/f noise
method = 4;

if method==1,
    % Growth rates for each species independent of spatial position
    r = comData.r_mean * abs(1+sqrt(comData.r_cv)*repmat(randn(comData.nSpecies,1),1,comData.nSites));
    % Rescale the growth rates
    r = r / mean(nonzeros(r));
elseif method==2,
    % For each resource define the position at which optimum resource occurs
    optimalX =  0.5*sqrt(comData.nResources)*randn(comData.nSpecies,comData.nResources);
    % weighting of each resource to the growth rate
    weights = 10*rand(comData.nSpecies,comData.nResources);


    % This lets growth rate be the product of N Gaussian distributions
    r = exp(- (sum(weights,2) * X.^2 ...
        - 2*sum(weights.*optimalX,2) * X + ...
        sum(weights.*optimalX.^2,2)*ones(size(X)))/comData.nResources);
    r_mean = mean(r,2);
    r = comData.r_mean * r ./ r_mean(:,ones(1,comData.nSites));
elseif method==3
    % Let the growth rate be a sum of linear contributions
    % First assumes that growth rate must be positive
    r_x0 = mean(max([0 comData.rMin]) + (comData.rMax - max([0 comData.rMin])) * rand(comData.nSpecies,comData.nResources),2);
    % Second assumes it can have the full range
    r_x1 = mean(comData.rMin + (comData.rMax - comData.rMin) * rand(comData.nSpecies,comData.nResources),2);

    % Place positive growth rate on left (X=1)
    resources.r_x0 = r_x0;
    resources.r_x1 = r_x1;

    ind = rand(comData.nSpecies,1)<0.5;
    % Half the time place positive growth rate on right (X=max(X))
    resources.r_x0(ind,:) = r_x1(ind,:);
    resources.r_x1(ind,:) = r_x0(ind,:);

    r = resources.r_x1(:,ones(size(R))) + ...
        ( resources.r_x0 - resources.r_x1 ) * (R - min(R))/(max(R)-min(R));
elseif method == 4,
    %    for i=1:comData.nResources,
    %        resources(i,:) = reshape(spatialPattern(comData.nSites,-4),1,prod(comData.nSites));
    %    end
    resources = 0;
%     r(1,:) = reshape(spatialPattern(comData.nSites,comData.spatialBeta),1,prod(comData.nSites));
%     r(1,:) = comData.rMin + (comData.rMax-comData.rMin) * (r(1,:)-min(r(1,:))) / (max(r(1,:)) - min(r(1,:)));
    for i=1:comData.nSpecies,
%       r(i,:) = r(1, :);
%       r(i,:) = r(1, :);
        r(i,:) = reshape(spatialPattern(comData.nSites,comData.spatialBeta),1,prod(comData.nSites));
        r(i,:) = comData.rMin + (comData.rMax-comData.rMin) * (r(i,:)-min(r(i,:))) / (max(r(i,:)) - min(r(i,:)));
    end
end

if simData.RandomSeed.growthRate~=0,
    rand('state',s);
    randn('state',sn);
end

return

%% Intra-DD gam
function [gam, intraDD] = intraDDGam(comData,simData,R),
% A function to give the itraspecific density dependence

% R is the resource coordinate
% Make sure it is a row vector
R = reshape(R,1,length(R));


if simData.RandomSeed.intraDD~=0,
    % Set random seed if required
    s = rand('state');
    rand('state',simData.RandomSeed.intraDD);
end


% Let the gam be a sum of linear contributions
% gam must be positive
gam_x0 = mean(comData.gamMin + (comData.gamMax - comData.gamMin) * rand(comData.nSpecies,comData.nResources),2);
gam_x1 = mean(comData.gamMin + (comData.gamMax - comData.gamMin) * rand(comData.nSpecies,comData.nResources),2);

intraDD.gam_x0 = gam_x0;
intraDD.gam_x1 = gam_x1;

gam = intraDD.gam_x1(:,ones(size(R))) + ...
    ( intraDD.gam_x0 - intraDD.gam_x1 ) * (R - min(R))/(max(R)-min(R));

if simData.RandomSeed.intraDD~=0,
    rand('state',s);
end


return

%% Interaction strengths
function I = interactionStrengths(webTopology,comData,simData),

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
function initialA = initialAbundance(comData,simData)

if simData.RandomSeed.initialAbundance~=0,
    % Set random seed if required
    s = rand('state');
    rand('state',simData.RandomSeed.initialAbundance);
    randn('state', simData.RandomSeed.initialAbundance);
end


initialA = comData.A_Init * ones(comData.nSpecies,prod(comData.nSites));
% Keep those species who have gone extinct
ind = initialA < comData.Amin;
initialA(ind) = 0;

% And remove species so that the community size on average is
% comData.nSpecies_mean
for k=1:prod(comData.nSites),
    if mod(k, round(comData.nSites/100)) == 0
      fprintf('Abundances done %g\n', k/prod(comData.nSites))
    end
    numSp = ceil(comData.nSpecies_mean*(1+randn(1)));
    if numSp>comData.nSpecies,
        numSp=comData.nSpecies;
    elseif numSp<0,
        numSp = 0;
    end
    rand_sp = randperm(comData.nSpecies);
    initialA(rand_sp(1:(comData.nSpecies-numSp)),k)= 0;
end

% % This just puts all species in the centre of the environment,
% % as a test of the spatial dispersal code
% initialA = zeros(size(initialA));
% ind = (comData.nSites(2)+1) * (round(comData.nSites(1)/2)-1)+1;
% initialA(:,ind) = 3*ones(comData.nSpecies,1);

% % This just puts all species at one random location,
% % as a test of the spatial dispersal code
% initialA = zeros(size(initialA));
% for i=1:comData.nSpecies,
%     ind = 1 + floor(rand*prod(comData.nSites));
%     initialA(i,ind) = 3;
% end

% Make A and c sparse matrices
initialA = sparse(initialA);

return
