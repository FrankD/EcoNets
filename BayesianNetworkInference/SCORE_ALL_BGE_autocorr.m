

function [SCORES] = SCORE_ALL_BGE_autocorr(data, T_0, T_m, v, alpha)

[n_nodes, n_obs] = size(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prior_0 = 1;
Prior_1 = 1/(n_nodes - 1);
Prior_2 = 1/((n_nodes -1)*(n_nodes -2)/2);
Prior_3 = 1/((n_nodes -1)*(n_nodes -2)*(n_nodes -3)/6);
Prior   = log([Prior_0, Prior_1, Prior_2, Prior_3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each network-node compute all local scores based on PRIOR:

n_nodes = n_nodes + 1;

for child = 1:(n_nodes-1)

% Score of the empty parent-set    
Counter = 1;    
Structure.parents(Counter,1:3) = [0 0 0];
log_likelihood  = Gauss_Score_complete_local(child, [n_nodes], ...
  data, T_0{1, child}, T_m{1, child}, v, alpha); %#ok<NBRAK>
Structure.scores(Counter) = log_likelihood  + Prior(1);
Counter = Counter +1;

% Local Prior scores of parent-sets of size 1:
for parent_1 = 1:(n_nodes - 1)
  if (parent_1 ~=child)
    Structure.parents(Counter,1:3) = [0, 0, parent_1];
    log_likelihood            = ...
      Gauss_Score_complete_local(child, [parent_1; n_nodes], ...
        data, T_0{1, child}, T_m{1, child}, v, alpha);
    Structure.scores(Counter) = log_likelihood  + Prior(2);
    Counter = Counter +1;
    
  end
end

% Local Prior scores of parent-sets of size 2:
for parent_1 = 1:(n_nodes-2)
  for parent_2 = (parent_1+1):(n_nodes-1)
    if (parent_1~=child && parent_2~=child)
      Structure.parents(Counter,1:3) = [0, parent_1, parent_2];
      log_likelihood            = ...
        Gauss_Score_complete_local(child,[parent_1;parent_2;n_nodes], ...
          data, T_0{1, child}, T_m{1, child}, v, alpha);
      Structure.scores(Counter) = log_likelihood + Prior(3) ;
      Counter = Counter +1;
      
    end
  end
end

% LocalPrior scores of parent-sets of size 3: 

for parent_1 = 1:(n_nodes-3)
  for parent_2 = (parent_1+1):(n_nodes-2)
    for parent_3 = (parent_2+1):(n_nodes-1)
      if (parent_1~=child && parent_2~=child && parent_3~=child)
        Structure.parents(Counter,1:3) = ...
          [parent_1, parent_2, parent_3];
          %[parent_1, parent_2, parent_3, n_nodes];
        log_likelihood            = ...
          Gauss_Score_complete_local(child, ...
            [parent_1;parent_2;parent_3;n_nodes], ...
            data, T_0{1, child}, T_m{1, child}, v, alpha);
        Structure.scores(Counter) = log_likelihood  + Prior(4);
        Counter = Counter +1;
        
      end
    end
  end
end

%highest = max(Structure.scores);

%Structure.scores = Structure.scores - highest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort with respect to the local Bayesian scores
[thrash, indicis] = sort(Structure.scores,'descend');

SCORES{child}.parents       = Structure.parents(indicis,:);
SCORES{child}.log_scores    = Structure.scores(indicis);
SCORES{child}.exp_scores    = exp(Structure.scores(indicis));

clear Structure;


end % "for child = 1:n_nodes"-loop

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [log_Score_local] = ...
  Gauss_Score_complete_local(index,index_parents, data, T_0, T_m, v, alpha)


[n, m] = size(data); 

relevant = sort([index;index_parents]);

[log_Score_i_nom]   = Gauss_Score_complete(length(relevant),      m, v, alpha, T_0(relevant,relevant)          , T_m(relevant,relevant));
     
[log_Score_i_denom] = Gauss_Score_complete(length(index_parents), m, v, alpha, T_0(index_parents,index_parents), T_m(index_parents,index_parents));
    
log_Score_local = log_Score_i_nom -log_Score_i_denom;     

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_Score_i] = Gauss_Score_complete(n_new, m, v, alpha, T_0, T_m)

if n_new==0
log_Score_i = 0;
else
      
sum_1 = (-1)*n_new*m/2 * log(2*pi) + (n_new/2) * log(v/(v+m));

sum_3 = (alpha/2) * log(det(T_0)) + (-1)*(alpha+m)/2 * log(det(T_m));

log_value_a = ( alpha    *n_new/2)*log(2) + (n_new*(n_new-1)/4)*log(pi);
log_value_b = ( (alpha+m)*n_new/2)*log(2) + (n_new*(n_new-1)/4)*log(pi);
 
for i=1:n_new
    log_value_a = log_value_a + gammaln(( alpha   +1-i)/2);
    log_value_b = log_value_b + gammaln(((alpha+m)+1-i)/2);
end
   
log_Score_i = sum_1 + (-1)*(log_value_a - log_value_b) + sum_3; 

end

return 
    









