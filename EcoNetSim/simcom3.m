% Script to look at the species community dynamics
% This version looks at a one dimensional environment
% This version also makes no changes to relative competitive ability

% Last edited 2/10/05  J. Yearsley

clear all
%close all

rand('state',sum(100*clock))
randn('state',sum(100*clock))

% Save model and results in different files
saveModel = 0;
saveResults = 0;
fileprefixModel = 'trial2Model';
fileprefixResults = 'trial2Results';

%% Parameter settings
% Simulation settings
% Set the random seeds for each random part of the model,
% zero means pick a truely random number
sim.RandomSeed.niche = 0; %4;
sim.RandomSeed.growthRate = 0; %[10,30,6];
sim.RandomSeed.intraDD = 0; % 9;
sim.RandomSeed.interactionStrengths = 0; % [13,22];
sim.RandomSeed.initialAbundance = 0;

sim.filename.model = 'notThere'; % File to load starting values for model if it exists
sim.filename.initialX = 'trialResults_18-Oct-2005#1'; % File to load starting initial abundance values

sim.repeat = 1; % Number of repeats
sim.iter = 1000;  % Number of iterations
sim.display = true;  % If set then display the results
sim.dt = 0.001; % Time step
sim.figureNum = 2;
sim.displayInterval = 100;
sim.movie = 0; % If true then create a movie of the simulation
sim.movieFile = [fileprefixResults '_movie_' date()];
sim.barPlot = 0; % If true then display bar chart, otherwise a line plot

% Community settings
com.nSites=100; % Number of spatial locations
com.nSpecies = 25; % Number of species
com.connectance = 0.1; % The connectance of the whole-web (assuming all species present)
com.Xmin = 0.01;  % Minimum abundance required for species survival

com.e_sigma = 0; % Variance of the environmental variation
com.interaction.PosMean = 5;
com.interaction.NegMean = 5;
com.interaction.PosSigma2 = 0.5;
com.interaction.NegSigma2 = 10;
com.gamMax = 10; % Max strength of intra-specific Gompertz density-dependence
com.gamMin = 10; % Min strength of intra-specific Gompertz density-dependence
%com.r_mean = 2; % Mean value of growth rate
com.rMax = 10; % Variation in growth rate
com.rMin = -10; % Variation in growth rate
com.rSigma = 0.1;
com.nResources = 1; % number of envrionmental variables that contribute to r
com.m = 0.1; % Mean distance of dispersal 

% Initial conditions
com.nSpecies_mean = 25;  % Mean number of species per site
com.X_Init = 1;

% Call script which will do the simulation
for run=1:sim.repeat,
    [X, deltaX, model, movieObj] = communitySimulation3(sim,com);

    if saveModel | saveResults,
        % Find a unique file name for saving data
        fileNum = 1;
        filenameModel = [fileprefixModel '_' date() '#' num2str(fileNum)];
        filenameResults = [fileprefixResults '_' date() '#' num2str(fileNum)];
        while exist([filenameModel '.mat'],'file') | exist([filenameResults '.mat'],'file'),
            fileNum = fileNum+1;
            filenameModel = [fileprefixModel '_' date() '#' num2str(fileNum)];
            filenameResults = [fileprefixResults '_' date() '#' num2str(fileNum)];
        end
    end
    if saveModel,
        save(filenameModel, 'model', 'com')
    end
    if saveResults,
        save(filenameResults, 'X','sim','com')
        %    if sim.movie,
        %      filename = [fileprefixResults '_movie_' date()];
        %      imwrite('figMovie',filename, 'bmp');
        %    end
    end
end
