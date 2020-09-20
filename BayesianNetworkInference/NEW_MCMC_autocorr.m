

function [Output] = NEW_MCMC_autocorr(step_iterations, dist, graph)

global SCORES;

global data;
global ns; 

n = length(SCORES);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Prior_0 = 1;
Prior_1 = 1/(n - 1);
Prior_2 = 1/((n -1)*(n -2)/2);
Prior_3 = 1/((n -1)*(n -2)*(n -3)/6);
Prior   = log([Prior_0, Prior_1, Prior_2, Prior_3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Counter = 1;
Output.dag{Counter} = graph; 


score_of_dag = 0;
    
for child=1:n
       
    parents_of_child = find(graph(:,child));
    log_likelihood = Gauss_Score_complete_local(child,parents_of_child);
    score_of_dag   = score_of_dag + Prior(length(parents_of_child)+1) + log_likelihood;

end
    
    
Output.score(Counter) = score_of_dag;




A = expm(graph') - eye(n);
A = (A>0);

%%%%%% Start the NEW-MCMC-ALGORITHM:

for t=1:step_iterations

  if(mod(t, round(step_iterations/100)) == 0)
    fprintf('Done %g%%', 100*(t/step_iterations));
  end

GRAPH_NEW = graph; % starting point of this MCMC-step
A_NEW     = A;

x = rand;  

        if (x<=(1/15))  % PERFORM A NEW REVERSAL MOVE 
            
 
        [x_edges, y_edges] = find(GRAPH_NEW);
        n_edges = sum(sum(GRAPH_NEW));
        
        
        if (n_edges==0)
            
           ACCEPTANCE = 1;
            
        else
        
        indicis = randperm(n_edges);
        index   = indicis(1);
        PARENT  = x_edges(index);
        CHILD   = y_edges(index);
        
        GRAPH{1}        = GRAPH_NEW; 
         
        GRAPH_NEW(:,PARENT)=0;  % PARENT and CHILD loose their parents
        GRAPH_NEW(:,CHILD)=0;
        
        A_NEW = expm(GRAPH_NEW') - eye(n);
        A_NEW = (A_NEW>0);
        
        LOG_SCORES_COPY = SCORES{PARENT}.log_scores;
        
        descendants_of_node = find(A_NEW(:,PARENT));
 
            for k = descendants_of_node'                                    
                % These nodes can NOT become parents of node!
                % EXCLUDE THEM:
                x1 = find(SCORES{PARENT}.parents(:,1)==k);
                x2 = find(SCORES{PARENT}.parents(:,2)==k);
                x3 = find(SCORES{PARENT}.parents(:,3)==k);
                x = [x1;x2;x3];
                LOG_SCORES_COPY(x) = -100000;
            end
            
            Indicator_Vector = -100000 * ones(length(LOG_SCORES_COPY),1);
                
               
            x = find(SCORES{PARENT}.parents(:,1)==CHILD | SCORES{PARENT}.parents(:,2)==CHILD | SCORES{PARENT}.parents(:,3)==CHILD);
            Indicator_Vector(x,1) = 0;  % indicates relevant parent-sets
               
                
            LOG_SCORES_COPY = LOG_SCORES_COPY + Indicator_Vector(:,1)'; % either + (-100000) or + 0
            
            
            C_NORM = max(LOG_SCORES_COPY);
            
            LOG_SCORES_COPY = LOG_SCORES_COPY - (C_NORM + 700);
            
            
            SCORES_COPY = exp(LOG_SCORES_COPY);
          
            cumulative = cumsum(SCORES_COPY)/sum(SCORES_COPY);
        
            x = rand;
            
           
            indicis = find((cumulative>=x));
     
            index = indicis(1); 
       
        new_parents_of_node = nonzeros(SCORES{PARENT}.parents(index,:));
        % Definitely including CHILD!!!
           
        GRAPH_NEW(new_parents_of_node,PARENT)=1;
  
        A_NEW = expm(GRAPH_NEW') - eye(n);
        A_NEW = (A_NEW>0);
            
            
        LOG_SCORES_COPY = SCORES{CHILD}.log_scores; 
         
        descendants_of_child = find(A_NEW(:,CHILD));
     
     
            for k = descendants_of_child'                           
            % These nodes can NOT become parents of child!
            % EXCLUDE THEM:
            x1 = find(SCORES{CHILD}.parents(:,1)==k);
            x2 = find(SCORES{CHILD}.parents(:,2)==k);
            x3 = find(SCORES{CHILD}.parents(:,3)==k);
            x = [x1;x2;x3];
            LOG_SCORES_COPY(x) = -100000;
            end
            
            C_NORM = max(LOG_SCORES_COPY);
            
            LOG_SCORES_COPY = LOG_SCORES_COPY - (C_NORM + 700);
            
            SCORES_COPY = exp(LOG_SCORES_COPY);
         
            cumulative = cumsum(SCORES_COPY)/sum(SCORES_COPY);
            
            x = rand;
   
            indicis = find((cumulative>=x));
       
            index = indicis(1); 
          
            
            new_parents_of_child = nonzeros(SCORES{CHILD}.parents(index,:));
        
            if ~isempty(new_parents_of_child)
            GRAPH_NEW(new_parents_of_child,CHILD)=1;
            end
        
            A_NEW = expm(GRAPH_NEW') - eye(n);
            A_NEW = (A_NEW>0);
        
            GRAPH{2}        = GRAPH_NEW; 
       
            m_parent = length(nonzeros(sum(GRAPH_NEW))); % number of nodes with parents
        
            ACCEPTANCE    = Acceptance_Probability_computation(PARENT,CHILD, GRAPH);
     
            ACCEPTANCE = min([ACCEPTANCE,1]);
            
            NEW_Indicator = 1;
        end
        
       
        
    else % PERFORM A STRUCTURE-MCMC-MOVE          
        
         
       
[neighbour, op, nodes, graph_neighbours] = neighbour_random(GRAPH_NEW, A_NEW, 3);

i = nodes(1); % the (modified) edge points from X_i  
j = nodes(2); % to X_j

% The local score of X_j (after operation):

       
        log_likelihood_new     = Gauss_Score_complete_local(j,find(neighbour(:,j))) + Prior(length(find(neighbour(:,j)))+1);
        LL_new                 = log_likelihood_new;

        log_likelihood_old     = Gauss_Score_complete_local(j,find(graph(:,j))) + Prior(length(find(graph(:,j)))+1);
        LL_old                 = log_likelihood_old;
            
        bf1 = exp(LL_new-LL_old); % Likelihood-Ratio(non-logarithmic) including non-uniform-prior

% If the operation is an edge-reversal...
if strcmp(op, 'rev')  
% ...also the local score of X_i changes: 

       
        log_likelihood_new     = Gauss_Score_complete_local(i,find(neighbour(:,i))) + Prior(length(find(neighbour(:,i)))+1);
        LL_new                 = log_likelihood_new;

      
        log_likelihood_old     = Gauss_Score_complete_local(i,find(graph(:,i)))+ Prior(length(find(graph(:,i)))+1) ;
        LL_old                 = log_likelihood_old;
            
        bf2 = exp(LL_new-LL_old); % Likelihood-Ratio(non-logarithmic) including non-uniform-prior

else
 bf2 = 1;
end

ACCEPTANCE = bf1 * bf2; % acceptance-probability
% without Hastings-Ratio -> Ratio of Bayesian Scores

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_neighbour = update_ancestor_matrix(A_NEW, op, i, j, neighbour);
[neighbour_neighbours] = n_neighbours(neighbour,A_neighbour, 3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Hasting_ratio = graph_neighbours/neighbour_neighbours;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACCEPTANCE = ACCEPTANCE * Hasting_ratio;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_NEW     = A_neighbour;
GRAPH_NEW = neighbour;

end       
    
 x = rand;
         
    if (x<=ACCEPTANCE)
        graph = GRAPH_NEW;
        A     = A_NEW;
        
    end
        
    
    if (mod(t,dist)==1)
    Counter = Counter+1;
    Output.dag{Counter}         = graph;
    end

end % loop: "for t=1:step_iterations"









return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTION

function ndx = subv2ind(siz, subv)

% siz can be a row or column vector of size d.
% subv should be a collection of N row vectors of size d.
% ind will be of size N * 1.
%
% Example:
% subv = [1 1 1;
%         2 1 1;
%         ...
%         2 2 2];
% subv2ind([2 2 2], subv) returns [1 2 ... 8]'
% i.e., the leftmost digit toggles fastest.

[ncases, ndims] = size(subv);

if all(siz==2)
  twos = pow2(0:ndims-1);
  ndx = ((subv-1) * twos(:)) + 1;
else
  cp = [1 ,cumprod(siz(1:end-1))']';
  ndx = (subv-1)*cp + 1;
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% update-ancestor-matrix-sub-functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = update_ancestor_matrix(A,  op, i, j, Nachbar)

switch op
case 'add',
 A = do_addition(A, op, i, j, Nachbar);
case 'del', 
 A = do_removal(A, op, i, j, Nachbar);
case 'rev', 
 A = do_removal(A, op, i, j, Nachbar);
 A = do_addition(A, op, j, i, Nachbar);
end



function A = do_addition(A, op, i, j, GRAPH_NEW)
% An edge from X_i to X_j has just been added 
% The ancestor-matrix A has to be updated

A(j,i) = 1; 
% Because X_i is an ancestor of X_j now

all_ancestors_i = find(A(i,:)); % All ancestors of X_i 

if ~isempty(all_ancestors_i) % If this set of ancestors of X_i is non-empty:  
 A(j,all_ancestors_i) = 1;   % These nodes become ancestors of X_j too 
end

all_ancestors_j   = find(A(j,:)); %  All ancestors of X_j (after adding X_i and its ancestors)
all_descendents_j = find(A(:,j)); %  All descendents of X_j (before adding X_i)

if ~isempty(all_ancestors_j)    % If the set of ancestors of X_j is non-empty:
    
    for k=all_descendents_j(:)' % For each descendent of X_j
    A(k,all_ancestors_j) = 1;   % add to the ancestors of this descendent k the ancestors of X_j 
    end
    
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = do_removal(A, op, i, j, GRAPH_NEW)
% An edge from X_i to X_j has just been deleted 
% The ancestor-matrix A has to be updated

% Determine all ancestors of X_j 
% and perform a topological sort of the nodes

descendents_of_j = find(A(:,j)); 
order            = topological_sort(GRAPH_NEW);

[THRASH, perm] = sort(order);

descendents_of_j_sorted = perm(descendents_of_j);

[THRASH, perm]   = sort(descendents_of_j_sorted);
descendents_of_j = descendents_of_j(perm);

% Update the descendents of X_j und fuer alle Nachfahren von Xj
A = update_row(A, j, GRAPH_NEW);

for k = descendents_of_j(:)'
A = update_row(A, k, GRAPH_NEW);
end

return;
%%%%%%%%%

function A = update_row(A, j, GRAPH_NEW)

% Determine the j-th row of A:
A(j, :) = 0; % Set this row = 0 
ps = find(GRAPH_NEW(:,j))'; % Determine the parents of X_j

if ~isempty(ps) % If X_j has parents
 A(j, ps) = 1; % All parents of X_j are andestors of X_j
end
for k=ps(:)' % For each parent k of X_j:
 ancestors_of_k = find(A(k,:)); % Determine the ancestors of k
 if ~isempty(ancestors_of_k)    % If there are some ancestors of k: 
   A(j, ancestors_of_k) = 1;    % set them ancestors of X_j
 end
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function order = topological_sort(dag) 

n          = length(dag);
indeg      = zeros(1,n);
zero_indeg = [];          % a stack of nodes with no parents

for i=1:n
indeg(i) = length(find(dag(:,i)));
  
     if indeg(i)==0
     zero_indeg = [i zero_indeg];
     end
  
end

t=1;
order = zeros(1,n);

while ~isempty(zero_indeg)
  v = zero_indeg(1); % pop v
  zero_indeg = zero_indeg(2:end);
  order(t) = v;
  t = t + 1;
  cs = find(dag(v,:));
  
  for j=1:length(cs)
    c = cs(j);
    indeg(c) = indeg(c) - 1;
      if indeg(c) == 0
      zero_indeg = [c zero_indeg]; % push c 
      end
  end
  
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ACCEPTANCE] = Acceptance_Probability_computation(node_1, node_2, GRAPH);

global SCORES;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GRAPH_1  = GRAPH{1};
GRAPH_2  = GRAPH{2};
n        = length(GRAPH_1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_edges_1  = sum(sum(GRAPH{1}));
n_edges_2  = sum(sum(GRAPH{2}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_parents_of_node_1  = find(GRAPH_2(:,node_1));   % is reached by the move
new_parents_of_node_2  = find(GRAPH_2(:,node_2));   % is reached by the move
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
node_a = node_2;
node_b = node_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_parents_of_node_a = find(GRAPH_1(:,node_a));    % for the inverse move
new_parents_of_node_b = find(GRAPH_1(:,node_b));    % for the inverse move
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+

LOG_SCORES_COPY_1 = SCORES{node_1}.log_scores;

GRAPH_1(:,node_1)=0;
GRAPH_1(:,node_2)=0;

A_1 = expm(GRAPH_1') - eye(n);
A_1 = (A_1>0);

descendants_1 =  find(A_1(:,node_1));
         
     for k = descendants_1'                                    
                % These nodes can NOT become parents of node!
                % EXCLUDE THEM:
                x1 = find(SCORES{node_1}.parents(:,1)==k);
                x2 = find(SCORES{node_1}.parents(:,2)==k);
                x3 = find(SCORES{node_1}.parents(:,3)==k);
                x = [x1;x2;x3];
                LOG_SCORES_COPY_1(x)=-100000;
     end
        
      
      Indicator_Vector = -100000 * ones(length(LOG_SCORES_COPY_1),1);
                              
      x = find(SCORES{node_1}.parents(:,1)==node_2 | SCORES{node_1}.parents(:,2)==node_2 | SCORES{node_1}.parents(:,3)==node_2);
      Indicator_Vector(x,1)=0;  % indicates relevant parent-sets
                      
      LOG_SCORES_COPY_1 = LOG_SCORES_COPY_1 + Indicator_Vector(:,1)';
      
 
      GRAPH_1(new_parents_of_node_1,node_1)=1;
         
       A_1 = expm(GRAPH_1') - eye(n);
       A_1 = (A_1>0);
        

      LOG_SCORES_COPY_2 = SCORES{node_2}.log_scores;
      
      descendants_2 = find(A_1(:,node_2));
      
 for k = descendants_2'                           
 % These nodes can NOT become parents of child!
 % EXCLUDE THEM:
 x1 = find(SCORES{node_2}.parents(:,1)==k);
 x2 = find(SCORES{node_2}.parents(:,2)==k);
 x3 = find(SCORES{node_2}.parents(:,3)==k);
 x = [x1;x2;x3];
 LOG_SCORES_COPY_2(x)= -100000;
 end
        
        GRAPH_2(:,node_a)   = 0;
        GRAPH_2(:,node_b)   = 0;
       
        A_2 = expm(GRAPH_2') - eye(n);
        A_2 = (A_2>0);
        
        descendants_a = find(A_2(:,node_a)); % all descendants of node_a
               
        LOG_SCORES_COPY_a = SCORES{node_a}.log_scores;
        
                for k=descendants_a'
                    x1 = find(SCORES{node_a}.parents(:,1)==k);
                    x2 = find(SCORES{node_a}.parents(:,2)==k);
                    x3 = find(SCORES{node_a}.parents(:,3)==k);
                    x = [x1;x2;x3];
                    LOG_SCORES_COPY_a(x) = -100000;
                end
    
            % Initialisation:
            Indicator_Vector = zeros(length(LOG_SCORES_COPY_a),1);
                   
            x = find(SCORES{node_a}.parents(:,1)==node_b | SCORES{node_a}.parents(:,2)==node_b | SCORES{node_a}.parents(:,3)==node_b );
            Indicator_Vector(x,1)=0;  % indicates relevant parent-sets
                   
            LOG_SCORES_COPY_a = LOG_SCORES_COPY_a + Indicator_Vector(:,1)';
           

             % Update: G_NEW and A_NEW 
       
             GRAPH_2(new_parents_of_node_a,node_a)=1;             
                
             A_2 = expm(GRAPH_2') - eye(n);
             A_2 = (A_2>0);           
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            descendants_b = find(A_2(:,node_b)); 
            % all descendants of node in G_2
            
            LOG_SCORES_COPY_b = SCORES{node_b}.log_scores;
          
                for k=descendants_b'
                    x1 = find(SCORES{node_b}.parents(:,1)==k);
                    x2 = find(SCORES{node_b}.parents(:,2)==k);
                    x3 = find(SCORES{node_b}.parents(:,3)==k);
                    x = [x1;x2;x3];
                    LOG_SCORES_COPY_b(x) = -100000;
                end
            
                
                
               
                 
            C_NORM_a = max(LOG_SCORES_COPY_a);
            C_NORM_b = max(LOG_SCORES_COPY_b);
            C_NORM_1 = max(LOG_SCORES_COPY_1);
            C_NORM_2 = max(LOG_SCORES_COPY_2);
            
            C_NORM = max([C_NORM_1,C_NORM_2,C_NORM_a,C_NORM_b]);      
                 
      LOG_SCORES_COPY_1 = LOG_SCORES_COPY_1 - (C_NORM + 700);            
      SCORES_COPY_1     = sum(exp(LOG_SCORES_COPY_1));
    
      LOG_SCORES_COPY_2 = LOG_SCORES_COPY_2 - (C_NORM + 700);                
      SCORES_COPY_2     = sum(exp(LOG_SCORES_COPY_2));
           
      LOG_SCORES_COPY_a = LOG_SCORES_COPY_a - (C_NORM + 700);      
      SCORES_COPY_a       = sum(exp(LOG_SCORES_COPY_a));
      

      LOG_SCORES_COPY_b = LOG_SCORES_COPY_b - (C_NORM + 700);
      SCORES_COPY_b     = sum(exp(LOG_SCORES_COPY_b));
      
      
      INDEX    = 0;
      INDEX_12 = 0;
      INDEX_ab = 0;
      
      if(SCORES_COPY_1 == 0 & SCORES_COPY_2==0)
          log_A_nominator   = 0;
          log_A_denominator = 100;
          INDEX=1;
      end
            
      if (SCORES_COPY_a == 0 & SCORES_COPY_b==0) 
          log_A_nominator   = 100;
          log_A_denominator = 0;
          INDEX=1;
      end
          
      if ((SCORES_COPY_1 == 0 | SCORES_COPY_2 == 0) & INDEX==0)  
           log_A_nominator = log(max([SCORES_COPY_1,SCORES_COPY_2])) + log(n_edges_1);
           INDEX_12=1;
      end    
         
      if ((SCORES_COPY_a == 0 | SCORES_COPY_b == 0) & INDEX==0)    
          log_A_denominator = log(max([SCORES_COPY_a,SCORES_COPY_b])) + log(n_edges_2);
          INDEX_ab=1;
      end 
          
      if (INDEX_12==0 & INDEX==0) 
          log_A_nominator = log(n_edges_1) + log(SCORES_COPY_1) + log(SCORES_COPY_2);
      end
      
      if(INDEX_ab==0 & INDEX==0)
          log_A_denominator = log(n_edges_2) + log(SCORES_COPY_a) + log(SCORES_COPY_b); 
      end
      
      
      
log_ACCEPTANCE = log_A_nominator - log_A_denominator;
ACCEPTANCE     = exp(log_ACCEPTANCE);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [graph_neighbours] = n_neighbours(G0, A, fan)

% Input: G0     - is the graph 
% Input: A      - is the ancestor-matrix of G0 
% Input: fan    - is the maximal fan-in 

% Output: The number of neighbour graphs

n_nodes = length(G0); % number of network-nodes

SS=sum(G0);   % Parent-set-size of each network-node (column-sums)
SS=(SS>=fan); % Nodes already having a parent-set of fan-in size (indicator-vector)
SS=find(SS);  % The indicis of these nodes 
             

% EDGE-DELETIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[I,J] = find(G0); % indicis of existing edges
E = length(I);    % number of exisiting edges

% EDGE-REVERSALS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = max(0, G0-(G0'*A)'); % L(i,j)=1 if and only if X_i->X_j is reversible 

L(SS,:)=0; % edges pointing from nodes with maximal fan-in are not reversible  
           
[IL, JL] = find(L);  % indices of reversible edges 

EL = length(IL); % number of reversible (non covered) 

% EDGE-ADDITIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gbar = ~G0;  % Alle non-existing edges 
% exclude: self-loops
self_loops = repmat(0, 1, n_nodes); 
diags = 1:n_nodes+1:n_nodes^2;
% Von 1 bis n^2 in (n+1)er-Schritten 
Gbar(diags) = self_loops;

GbarL = Gbar-A; % avoid cycles 
GbarL(:,SS)=0;  % exclude edges leading to nodes with maximal fan-in      
[IbarL, JbarL] = find(GbarL);  % indicis of edges that can be added 
EbarL = length(IbarL); % number of these edges

graph_neighbours = E+EL+EbarL;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Nachbar, op, nodes, graph_neighbours] = neighbour_random(G0, A, fan)

% Input: G0     - is the graph 
% Input: A      - is the ancestor-matrix of G0 
% Input: fan    - is the maximal fan-in 

% Output: A randomly choosen neighbour-graph of G0

n_nodes = length(G0); % number of network-nodes

SS=sum(G0); % parent-set-size of each network-node (column-sums)

SS=(SS>=fan); % Nodes already having a parent-set of fan-in size (indicator-vector)
SS=find(SS); % The indicis of these nodes 
             

% EDGE-DELETIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[I,J] = find(G0); % indicis of existing edges
E = length(I);    % number of exisiting edges


   

% EDGE-REVERSALS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = max(0, G0-(G0'*A)'); % L(i,j)=1 if and only if X_i->X_j is reversible 

L(SS,:)=0; % edges pointing from nodes with maximal fan-in are not reversible  
           
[IL, JL] = find(L);  % indices of reversible edges 

EL = length(IL); % number of reversible (non covered) 


% EDGE-ADDITIONS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gbar = ~G0;  % Alle non-existing edges 
% exclude: self-loops
self_loops = repmat(0, 1, n_nodes); 
diags = 1:n_nodes+1:n_nodes^2;
% Von 1 bis n^2 in (n+1)er-Schritten 
Gbar(diags) = self_loops;


GbarL = Gbar-A; % avoid cycles 
GbarL(:,SS)=0;  % exclude edges leading to nodes with maximal fan-in
            
[IbarL, JbarL] = find(GbarL);  % indicis of edges that can be added 
EbarL = length(IbarL); % number of these edges

graph_neighbours = (E+EL+EbarL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Randomly chose a neighbour graph:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bestimme zunaechst, die Operation: del, rev oder add
Probs=[E, EL, EbarL]/(graph_neighbours);

R=rand(1); % Eine (0,1)-Zufallszahl, gleichverteilt

% Bestimme die Verteilungsfunktion von Probs:
cumprob = cumsum(Probs(:));
cumprob = cumprob(1:end-1);

% Wieviele Werte der Verteilungsfunktion werden ueberschritten?
Zufall = sum(R > cumprob)+1;
      
% Beruecksichtigt wurde, dass durch die drei Operationen 
% jeweils unterschiedlich viele Nachbar-Graphen erreicht werden koennen!
% Die Variable "Zufall" ist aus der Menge {1,2,3}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schreibe die Adjacency-Matrix als Spaltenvektor:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gkopie = G0(:);

% In Abhaengigkeit von der gewaehlten Operation:

if Zufall==1 % EDGE-DELETION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
op='del';
%Auswahl=random('Discrete Uniform',E,1,1);
Auswahl=floor(rand()*E) + 1;
% Waehle zufaellig eine der E Kanten, die entfernt werden kann
x = I(Auswahl); % Die Indices
y = J(Auswahl); % dieser Kante

ndx = subv2ind([n_nodes n_nodes], [x y]); % Berechne den 1-dim Index dieser Kante
% Fasse dazu [x y] als Element einer nxn-Matrix auf

Gkopie(ndx)=0; % Enferne die Kante
Nachbar = reshape(Gkopie, [n_nodes n_nodes]); % Forme die neue Adjacency-Matrix

elseif Zufall==2 % EDGE-REVERSAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
op='rev';
%Auswahl=random('Discrete Uniform',EL,1,1);
Auswahl=floor(rand()*EL) + 1;
x = IL(Auswahl); % Indices der Kante
y = JL(Auswahl); % die umgekehrt wird

ndx = subv2ind([n_nodes n_nodes], [x y]); 
rev_ndx = subv2ind([n_nodes n_nodes], [y x]); % die 1-dim Indices
Gkopie(ndx) = 0; % Kante entfernen
Gkopie(rev_ndx) = 1; % neue Kante hinzufuegen
Nachbar = reshape(Gkopie, [n_nodes n_nodes]); % wieder eine Adjacency-Matrix formen

else % EDGE-ADDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
op='add';
%Auswahl=random('Discrete Uniform',EbarL,1,1);
Auswahl=floor(rand()*EbarL) + 1;
x = IbarL(Auswahl);
y = JbarL(Auswahl);
 
ndx = subv2ind([n_nodes n_nodes], [x y]);
Gkopie(ndx) = 1;
Nachbar = reshape(Gkopie, [n_nodes n_nodes]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodes=[x y];


return;
