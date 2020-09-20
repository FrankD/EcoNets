

function [Scores] = COMPUTE_SCORES_FROM_Output_autocorr(Output)

n_sample = length(Output.dag);


n_nodes = length(Output.dag{1});

Scores = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prior_0 = 1;
Prior_1 = 1/(n_nodes - 1);
Prior_2 = 1/((n_nodes -1)*(n_nodes -2)/2);
Prior_3 = 1/((n_nodes -1)*(n_nodes -2)*(n_nodes -3)/6);
Prior   = log([Prior_0, Prior_1, Prior_2, Prior_3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n_sample
    
    graph = Output.dag{i};
        
    Score_i = 0;
    
    for j=1:n_nodes
        Score_i = Score_i + Gauss_Score_complete_local(j,find(graph(:,j))) + Prior(length(find(graph(:,j)))+1);
    end
    
    Scores = [Scores,Score_i];
    
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [log_Score_local] = Gauss_Score_complete_local(index,index_parents)

global data;
global T_0;
global T_m;
global v;
global alpha;

T_0_index = T_0{1, index};
T_m_index = T_m{1, index};

[n, m] = size(data); 

n = n +1;

index_parents = [index_parents; n];

relevant = sort([index;index_parents]);

[log_Score_i_nom]   = Gauss_Score_complete(length(relevant),      m, v, alpha, T_0_index(relevant,relevant)          , T_m_index(relevant,relevant));
     
[log_Score_i_denom] = Gauss_Score_complete(length(index_parents), m, v, alpha, T_0_index(index_parents,index_parents), T_m_index(index_parents,index_parents));
    
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
