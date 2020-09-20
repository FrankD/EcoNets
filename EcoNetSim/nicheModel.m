function [P,trophicLevel] = nicheModel(S,C,redundancy,randState)
% A function which generates a food web following the niche model
% of Williams & Martinez (2000)"Simple rules yield complex food webs"
% Nature, Vol 404, p 180-183
%
% S  the number of species in the web
% C  the connectance of the web
% redundancy  (optional argument) if true then don't remove trophically identical species
%
% The script generates a prey matrix specifying the prey of each species
% P(:,3) gives the prey of the third species
% There is at least one basal species.
% trophicLevel is an optional output which gives the trophic level!!

% Written by Jon Yearsley 15/04/2005  j.yearsley@macaulay.ac.uk

%% Check input arguments
if nargin<3,
    redundancy = false;
end

disp('Start Generating Niche Model')

if C<0 || C>0.5
    error('Connectance must be greater than zero and less than 0.5')
end

if S<1
    error('The number of species in a web must be positive')
end

if round(S)-S>eps
    S = round(S);
    warning('The number of species is not interger. Rounding to nearest integer.')
end

if S==1,
    P=sparse(1,1,0);
    trophicLevel = 1;
    return;
end

if nargin > 3 & randState~=0,
    % Set random seed if required
    s = rand('state');
    rand('state',randState);
end

%% Start to simulate a web
checkSpecies = true;
continuous = false;
disp('Iterate while not continuous');
i = 0;
while ~continuous || any(checkSpecies)

    i = i + 1;
    
    if ~continuous,
        % Generate parameters for the model
        checkSpecies = true(S,1);
        [n, niche] = paramGenerate(S,C);
    end

    % Generate predation matrix
    P = webGenerate(n,niche,S,C);

    % Check to see if web is continuous (i.e. not disjoint)
    continuous = ( nnz((P+P')^S) == S^2 );

    % Number of prey for each species
    prey = full(sum(P,1));

    if continuous && ~redundancy,
        % Check to see if two species are trophically identical
        checkSpecies(prey==0) = false;  % ignore basal species
        
        if any(checkSpecies)
            % Just look at species which weren't independent before
            % Sort these species in terms of the number of prey
            preySort = sort(prey(checkSpecies));
            preyItems = [preySort(1) preySort(find(diff(preySort))+1)];
            for p=1:length(preyItems),
                ind = find(prey==preyItems(p));
                % Find species which are trophically independent
                indep = all(abs(corrcoef(P(:,ind))-speye(length(ind))-1)>eps*10);
                % No longer consider species that are trophically independent
                checkSpecies(ind(indep)) = false;
                % Finally if the dependent species were all previously
                % dependent, ignore one, and resimulate the others
                if all(checkSpecies(ind(~indep))),
                    % No longer consider the first of these dependent species
                    checkSpecies(cumsum(checkSpecies(ind(~indep)))==1) = 0;
                end
            end
            % Resimulate the model parameters for the dependent species
            [n_new,niche_new] = paramGenerate(nnz(checkSpecies),C);
            n(checkSpecies) = n_new;
            niche(checkSpecies,:) = niche_new;
        end
    end
    if redundancy,
        checkSpecies = false;
    end
end

%% Finally sort web in order of trophic level
% Minimum trophic level is the minimum food chain length below the species
trophicLevel = -ones(S,1);
level = 0;
trophicLevel(prey==0,:) = level; % These are basal species

% Next line gives minimum food chain length below a species
while any(trophicLevel<0),
    level = level+1;
    trophicLevel(any(P(trophicLevel==level-1,:)) & trophicLevel'<0) = level;
end

% % Next line gives maximum food chain length below a species
% while any(trophicLevel<0),
%     level = level+1;
%     % If all the prey of a species already has a tropic level
%     trophicLevel(all(P(trophicLevel~=-1,:),1)) = level;
% end

[trophicLevel,ind] = sort(trophicLevel);
P = P(ind,ind);

if nargin==4 & randState~=0,
    % Reset random number seed to its state before the function call
    rand('state',s)
end

return

%%%%%%%%%%%%%%%%
% Below are sub functions used in the script
%%%%%%%%%%%%%%%%
%% Generate niche model parameters
function [n,niche] = paramGenerate(S,C)

% Generate a uniform random number for each species and sort in ascending
% order
n = sort(rand(S,1));

% Generate a beta distribution random deviate (with ALPHA=1) for each
% species with expectation=2C
% r is the niche width of a species
BETA = 1/(2*C)-1;
r = 1 - (1-rand(S,1)).^(1/BETA);

% Make the species with the smallest n have zero r, to ensure one basal
% species
r(1) = 0;

% Generate a uniform random deviate from r/2 to n for each species
% This is the central position of the niche
c = r/2 + (n-r/2).*rand(S,1);

% Calculate niche position
niche = [n - c - r/2 n - c + r/2];

return

%%%%%%%%%%%%%%%%
%% Function that generates web from niche model parameters
function P = webGenerate(n,niche,S,C)

P = spalloc(S,S,ceil(S^2*C));
for s=1:S,
    P(:,s) = ( n>niche(s,1) & n<niche(s,2) );
end
