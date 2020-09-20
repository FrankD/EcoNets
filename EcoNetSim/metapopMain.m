% Script to enter stuff into metapop script

% Generate resources
clear all

nSites = 1;
nSpecies = 5;
nResources = 2;
nTermsR = 2;
nTermsP = 2;

degreeR = 2;
degreeP = 2;

if nSites == 1,
    R = randn(1,nResources);
else
    R = zeros(nSites,nResources);
    for i=1:nResources,
        R(:,i) = spatialPattern([nSites,1],-2);
    end
end

betaR = zeros(nResources,nTermsR,nSpecies);
betaP = zeros(nSpecies,nTermsP+1,nSpecies);
for i=1:nSpecies,
    for term=1:nTermsR,
        % Select a random indeces for each term
        ind = randperm(nResources);
        betaR(ind(1:degreeR),term,i) = 1;
    end

    for term=1:nTermsP,
        % Select a random indeces for each term
        ind = randperm(nSpecies);
        betaP(ind(1:degreeP),term,i) = 1;
    end
    betaP(i,nTermsP+1,i) = 1; % Self interaction
end

wR = ones(nSpecies,nTermsR);
wP = ones(nSpecies,nTermsP+1);
wP(:,end) = -1;

x0 = zeros(nSites,nSpecies);
[xSol,fval,exitflag] = fsolve(@(x) metapop(x,R,betaR,betaP,wR,wP),x0);

pSol = exp(xSol) / (1+exp(xSol));
