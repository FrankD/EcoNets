% Function bayes_edge_autocorr uses MCMC with a special edge reversal move
% to determine the posterior edge probabilites of a Bayesian network from
% the data. The MCMC algorithm is the on described in,
%
% Grzegorczyk, M.; and Husmeier, D. (2008): Improving the structure MCMC 
% sampler for Bayesian networks by introducing a new edge reversal move.
% Machine Learning, 71, 265-305. 
% 
% This algorithm is further extended to produce a network that can deal
% with data where spatial autocorrelation (autocorrelation between adjacent
% or close data points) is a significant problem.
% 
% This is achieved conceptually by considering each node in the network to 
% be linked to a node that represents the autocorrelation. Thus, the effect
% from autocorrelation is accounted for, and the MCMC algorithm can
% concentrate on other effects.
%
% From a practical point of view, it was easiest to simply extend the
% precomputed scores so that there are N sets of scores (one for each
% original node), each with n+1 nodes, where the last node of score set i
% represents the autocorrelation for node i. In the MCMC algorithm, the
% right score set is chosen depending on which node we are currently
% modifying the parentset of.
%
% Inputs:
% local_data - The data matrix, an unnormalised NxM matrix, where N is the
% number of nodes and M is the number of locations/data points
% iterations - The total number of iterations for the MCMC algorithm
%
% Outputs:
% uge_scores - Undirected posterior edge probabilities 
% dge_scores - Directed posterior edge probabilites
% DAGs - The networks that were sampled during MCMC (one ever 100
% iterations)
% Scores - Likelihoods of the networks that were sampled
%
% Author: Frank Dondelinger (frank.dondelinger@gmail.com)
% Last Modified: November 2008
function [uge_scores, dge_scores, DAGs, Scores] = ... 
  bayes_edge_autocorr(local_data, iterations)

% Normalise Data
global data;
data = standardizeCols(local_data')';

n_nodes = size(data, 1);

global T_0;
global T_m;
global v;
global alpha;

% Generate Prior info
T_0 = cell(1, n_nodes);
T_m = cell(1, n_nodes);

for node = 1:n_nodes
  [T_0_node, T_m_node, v, alpha] = ...
    Compute_Prior_Info_BGE_autocorr(data, node);
  
  T_0{1, node} = T_0_node;
  T_m{1, node} = T_m_node;
end

% Calculate Scores
global SCORES;
[SCORES] = SCORE_ALL_BGE_autocorr(data, T_0, T_m, v, alpha);

% Initial Graph
graph = zeros(n_nodes, n_nodes);

step_iterations = iterations;
dist = round(iterations/100);

% Do MCMC
[Output] = NEW_MCMC_autocorr(step_iterations, dist, graph);

DAGs = Output;

display 'Scores for DAGs';

% 6. The Output variable does not contain the scores of the DAGs
% Compute the scores:
[Scores] = COMPUTE_SCORES_FROM_Output_autocorr(Output); 

display 'Edge Scores';
% Compute Posterior probabilities for edges
[uge_scores, dge_scores] = AUROC(Output, round(size(Output.dag, 2)/10));
