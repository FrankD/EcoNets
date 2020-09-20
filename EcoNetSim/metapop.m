function dlogPdt = metapop(x,R,betaR,betaP,wR,wP),

P = exp(x) ./ (1+exp(x));

[nSites,nSpecies] = size(x);
[nResources,nTermsR,nSpecies] = size(betaR);
[nSpecies,nTermsP,nSpecies] = size(betaP);

for i=1:nSpecies
    % The squeeze removes the spatial dimension if there's only one site
    R_i = squeeze(repmat(R,[1,1,nTermsR]));
    P_i = squeeze(repmat(P,[1,1,nTermsP]));
    
    if nSites==1,
        
        betaR_i = betaR(:,:,i);
        betaP_i = betaP(:,:,i);
        
        dlogPdt(:,i) = sum(wR(i,:) .* prod(R_i.^betaR_i , 1))/sum(wR(i,:));
        
        dlogPdt(:,i) = dlogPdt(:,i) + ...
            sum(wP(i,:) .* prod(P_i.^betaP_i ,1) ) / sum(wP(i,:));
        
    else
        wR_i = repmat(wP(i,:),[nSites,1]);
        wP_i = repmat(wP(i,:),[nSites,1]);
        
        betaR_i = permute(repmat(betaR(:,:,i),[1,1,nSites]),[3 1 2]);
        betaP_i = permute(repmat(betaP(:,:,i),[1,1,nSites]),[3 1 2]);
        
        dlogPdt(:,i) = squeeze(sum(wR_i .* squeeze(prod(R_i.^betaR_i , 2)),2))/sum(wR(i,:));
        
        dlogPdt(:,i) = dlogPdt(:,i) + ...
            squeeze(sum(wP_i .* squeeze(prod(P_i.^betaP_i , 2)) , 2))/sum(wP(i,:));
    end
end
