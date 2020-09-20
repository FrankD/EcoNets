% Script to look at the species community dynamics
% This version looks at a one dimensional environment
% This version also makes no changes to relative competitive ability

% Last edited 2/10/05  J. Yearsley

clear all
%close all

rand('state',sum(100*clock))
randn('state',sum(100*clock))

%P = path();
%path(P,'/data1/jl5362/JonYStuff')

% Save model and results in different files
saveModel = 1;
saveResults = 1;
fileprefix = 'simulatedDataToy';
%fileprefixDispersalMat = 'trial11x11Model';
%fileprefixResults = 'trial11x11Results';

%% Parameter settings
% Simulation settings
% Set the random seeds for each random part of the model,
% zero means pick a truely random number
sim.RandomSeed.niche = 0; %10;
sim.RandomSeed.growthRate = 0; %[10,30,6];
sim.RandomSeed.intraDD = 0; %5;
sim.RandomSeed.interactionStrengths = 0; %[13,22];
sim.RandomSeed.initialAbundance = 0;

sim.filename.model = 'notThere'; % File to load starting values for model if it exists
sim.filename.DispersalMatrix = 'dispersalMatrix10x10'; % File to load starting values for model if it exists
sim.filename.initialA = 'trial31x31_Results_25-Oct-2005#1'; % File to load starting initial abundance values

sim.repeat = 10; % Number of repeats
sim.iter = 100;  % Number of iterations
sim.display = 0;  % If set then display the results
sim.debug = 0;
sim.dt = 0.001; % Time step
sim.figureNum = 2;
sim.displayInterval = 1;
sim.movie = 0; % If true then create a movie of the simulation
sim.movieFile = [fileprefix '_movie_' date()];
sim.displayAmax = 3;
% displayCode = 0   Plot species richness
% displayCode = -1   Plot mean abundance
% displayCode = -2  Plot identity of the most abundant
% displayCode = -4  Plot identity of the most abundant
% displayCode = -3  Plot an abundance transect through the data
% displayCode = i   Plot abundance of species i
sim.displayCode = 10;

sim.changeGrowthRate = 0;
sim.changeInteractions = 0;

% Community settings
com.nSites=[5,5]; % Number of spatial locations
com.nSpecies = 3; % Number of species
com.connectance = 0.1; % The connectance of the whole-web (assuming all species present)
com.Amin = 0.05;  % Expected abundance required for species survival (exponential distribution)

com.e_sigma = 0; % Variance of the environmental variation
com.interaction.PosMean = 5;
com.interaction.NegMean = 5;
com.interaction.PosSigma2 = 0.5;
com.interaction.NegSigma2 = 0.5;
com.gamMax = 10; % Max strength of intra-specific Gompertz density-dependence
com.gamMin = 10; % Min strength of intra-specific Gompertz density-dependence
%com.r_mean = 2; % Mean value of growth rate
com.rMax = 10; % Variation in growth rate
com.rMin = -10; % Variation in growth rate
com.nResources = 1; % number of envrionmental variables that contribute to r
com.spatialBeta = -4;
com.m = 0.1; % Mean distance of dispersal

% Initial conditions
com.nSpecies_mean = 3;  % Mean number of species per site
com.A_Init = 3;

% Call script which will do the simulation
for run=1:sim.repeat,
    [A, model, movieObj] = communitySimulation_2d_ts(sim,com);

    % Finally calculate some things
    %richness = full(reshape(sum(A>com.Amin,1),com.nSites));
%    tmp = quantile(A,[0,0.25,0.5,0.75,1],1);
%    Amin = reshape(tmp(1,:),com.nSites);
%    Alquart = reshape(tmp(2,:),com.nSites);
%    Amed = reshape(tmp(3,:),com.nSites);
%    Auquart = reshape(tmp(4,:),com.nSites);
%    Amax = reshape(tmp(5,:),com.nSites);
%    [dum,tmp] = max(A,[],1);
%    domID = reshape(tmp,com.nSites);


    if saveModel | saveResults,
        % Find a unique file name for saving data
        fileNum = 1;
        filenameModel = [fileprefix '_Model_' date() '#' num2str(fileNum)];
        filenameResults = [fileprefix '_Results_' date() '#' num2str(fileNum)];
        filenameDispersalMatrix = [fileprefix '_DispersalMatrix_' date() '#' num2str(fileNum)];
        while exist([filenameModel '.mat'],'file') | ...
                exist([filenameResults '.mat'],'file') | ...
                exist([filenameDispersalMatrix '.mat'],'file'),

            fileNum = fileNum+1;
            filenameModel = [fileprefix '_Model_' date() '#' num2str(fileNum)];
            filenameResults = [fileprefix '_Results_' date() '#' num2str(fileNum)];
            filenameDispersalMatrix = [fileprefix '_DispersalMatrix_' date() '#' num2str(fileNum)];
        end
    end
    if saveModel,
        save(filenameModel, 'model', 'com')
    end
    if saveResults,
        save(filenameResults, 'A','sim','com','richness')

        dispersalMatrix = model.dispersalMatrix;
        xCoord = model.xCoord;
        yCoord = model.yCoord;
        save(filenameDispersalMatrix, 'dispersalMatrix','xCoord','yCoord')

        if sim.movie,
            filenameMovie = [fileprefix '_MovieData_' date() '#' num2str(fileNum)];
            save(filenameMovie,'movieObj','com','sim')
        end
    end
end


