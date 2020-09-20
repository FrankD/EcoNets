% Script to look at the species community dynamics
% This version looks at a two dimensional environment
% This version also makes no changes to relative competitive ability

% Last edited 2/10/05  J. Yearsley

function [model, A, extinct] = generateSimData(num_species, loc_dim, steps, ...
  connectance, pos_interactions, neg_interactions, gamma, noise, rate, ...
  beta)

%close all

rand('state',sum(100*clock))
randn('state',sum(100*clock))

%P = path();
%path(P,'/data1/jl5362/JonYStuff')


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

sim.repeat = 1; % Number of repeats
sim.iter = steps;  % Number of iterations
sim.display = 0;  % If set then display the results
sim.debug = 0;
sim.dt = 0.01; % Time step
sim.figureNum = 2;
sim.displayInterval = 1;
sim.movie = 0; % If true then create a movie of the simulation
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
com.nSites= loc_dim; % Number of spatial locations
com.nSpecies = num_species; % Number of species
if(nargin < 4)
  com.connectance = 0.1; % The connectance of the whole-web (assuming all species present)
else
  com.connectance = connectance;
end
com.Amin = 0.05;  % Expected abundance required for species survival (exponential distribution)


if(nargin < 8)
  com.e_sigma = 0; % Variance of the environmental variation
else
  com.e_sigma = noise;
end
% Good parameters 0, 10 with gam 8, growth rate 10, 0, connectance 0.1
if(nargin < 5)
  com.interaction.PosMean = 4;
else
  com.interaction.PosMean = pos_interactions;
end

if(nargin < 6)
  com.interaction.NegMean = 8;
else
  com.interaction.NegMean = neg_interactions;
end

com.interaction.PosSigma2 = 0.5;
com.interaction.NegSigma2 = 0.5;
if(nargin < 7)
  com.gamMax = 13;
  com.gamMin = 13;
else
  com.gamMax = gamma; % Max strength of intra-specific Gompertz density-dependence
  com.gamMin = gamma; % Min strength of intra-specific Gompertz density-dependence
end
%com.r_mean = 2; % Mean value of growth rate
if(nargin < 9)
  com.rMax = 10; % Variation in growth rate
else
  com.rMax = rate; % Variation in growth rate
end

com.rMin = 0; % Variation in growth rate
com.nResources = 1; % number of envrionmental variables that contribute to r

if(nargin < 10)
  com.spatialBeta = -1;
else
  com.spatialBeta = beta;
end

com.m = 0.1; % Mean distance of dispersal

% Initial conditions
com.nSpecies_mean = num_species/2;  % Mean number of species per site
com.A_Init = 3;

% Call script which will do the simulation
[A, model, movieObj] = communitySimulation_2d(sim,com);

[A, model, extinct] = remove_extinct(full(A), model);
%extinct = 0;

end


function [A, model, extinct] = remove_extinct(A_full, model)
%% Function to find and remove extinct species

extinct_inds = abs(sum(A_full,2))<eps;
extinct = find(extinct_inds);

A = A_full(~extinct_inds,:);

end