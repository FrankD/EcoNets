
function [UGE, DGE] = AUROC(Run, burn)

DAGS = Run.dag;

% The number of sampled DAGs:
n_samples = length(DAGS);

% The number of network-nodes
n_nodes = size(DAGS{1},1);

% Initialisation for relation-feature aposteriori-probabilities:
DIRECTED   = zeros(n_nodes,n_nodes);
UNDIRECTED = zeros(n_nodes,n_nodes);

for i = (burn+2):n_samples
        
[Current_DAG]   = DAGS{i};

[Current_CPDAG] = dag_to_cpdag(Current_DAG);

Mat_fix = zeros(n_nodes,n_nodes); % count compelled and non-
Mat_rev = zeros(n_nodes,n_nodes); % compelled edges in these matrices

indicis_fix = find(Current_CPDAG>0);
indicis_rev = find(Current_CPDAG<0);

Mat_fix(indicis_fix)=1;
Mat_rev(indicis_rev)=1;

DIRECTED_new = Mat_fix + (Mat_rev + Mat_rev');

DIRECTED = DIRECTED + DIRECTED_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And now determine: (Undirected) Edges-Relation-Features:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UNDIRECTED_new = (Current_DAG + Current_DAG');

UNDIRECTED     =  UNDIRECTED + UNDIRECTED_new;

end % end of loop "for i=1:n_samples"

n_samples = n_samples - (burn+1);

% OUTPUT:
DGE  = DIRECTED/n_samples;
UGE = UNDIRECTED/n_samples;


return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [label] = dag_to_cpdag(dag)
% LABEL-EDGES-Algorithms
% 	+1 if the edge is "compelled" 
%	-1 if the edge is "reversible"

N=length(dag);
[order xedge yedge] = order_edges(dag); 

label = 2*dag; 

NbEdges = length(xedge) ; 


for Edge=1:NbEdges 
    xlow=xedge(Edge); 
    ylow=yedge(Edge); 
    
    if label(xlow,ylow)==2 
        fin = 0; 
        wcompelled = find(label(:,xlow)==1); 
        parenty = find(label(:,ylow)~=0); 
                                         

        for w = wcompelled' 
           
            if ~ismember(w,parenty) 
                label(parenty,ylow)=1; 
                fin = 1;  
            elseif fin == 0 
                
                label(w,ylow)=1; 
               
            end
            
        end
        
        if fin == 0
            parentx = [xlow ; find(label(:,xlow)~=0)];
            
            if ~isempty(mysetdiff(parenty,parentx)) 
                label(find(label(:,ylow)==2),ylow)=1;  
            else	
                label(xlow,ylow)=-1; % xlow->ylow "Reversible" 
                %label(ylow,xlow)=-1; 
                label(find(label(:,ylow)==2),ylow)=-1; 
                
            end
        end
    end
end

%%%========================================================================================
function [order, x, y] = order_edges(dag)

N=length(dag);  
order = zeros(N,N); 

node_order = topological_sort(dag); 
[tmp oo] = sort(node_order);

dag=dag(node_order,node_order); 

[x y]=find(flipud(dag)==1);

nb_edges=length(x); 

if nb_edges~=0  
  order(sub2ind([N N],N+1-x,y))=1:nb_edges ; 
end

order=order(oo,oo); 

x=node_order(N+1-x);
y=node_order(y);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function order = topological_sort(A)
% TOPOLOGICAL_SORT Return the nodes in topological order (parents before children).
% order = topological_sort(adj_mat)

n = length(A);
indeg = zeros(1,n);
zero_indeg = []; % a stack of nodes with no parents
for i=1:n
  indeg(i) = length(parents(A,i));
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
  cs = children(A, v);
  for j=1:length(cs)
    c = cs(j);
    indeg(c) = indeg(c) - 1;
    if indeg(c) == 0
      zero_indeg = [c zero_indeg]; % push c 
    end
  end
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ps = parents(adj_mat, i)
% PARENTS Return the list of parents of node i
% ps = parents(adj_mat, i)

ps = find(adj_mat(:,i))';

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = mysetdiff(A,B)
% MYSETDIFF Set difference of two sets of positive integers (much faster than built-in setdiff)
% C = mysetdiff(A,B)
% C = A \ B = { things in A that are not in B }

if isempty(A)
  ma = 0;
else
  ma = max(A);
end

if isempty(B)
  mb = 0;
else
  mb = max(B);
end

if ma==0 
  C = [];
elseif mb==0
  C = A;
else % both non-empty
  %bits = sparse(1, max(ma,mb));
  bits = zeros(1, max(ma,mb));
  bits(A) = 1;
  bits(B) = 0;
  C = A(logical(bits(A)));
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cs = children(adj_mat, i, t)
% CHILDREN Return the indices of a node's children in sorted order
% c = children(adj_mat, i, t)

if nargin < 3 
  cs = find(adj_mat(i,:));
else
  if t==1
    cs = find(adj_mat(i,:));
  else
    ss = length(adj_mat)/2;
    j = i+ss;
    cs = find(adj_mat(j,:)) + (t-2)*ss;
  end
end

return



 