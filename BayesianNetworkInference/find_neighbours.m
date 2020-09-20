% Function all_neighbours computes the direct neighbours of each location,
% assuming there are S*S locations in total, arranged on an SxS grid.
% 
% Inputs
% data - The original data, an NxM matrix, where N is the number of species
% and M is the number of locations.
%
% Outputs:
% all_neighbours - A 1xM cell array containing the neighbours for each
% location. all_neighbours{i} is a matrix of size 1xL, where L is the
% number of neighbours of location i.
%
% Author: Frank Dondelinger (frank.dondelinger@gmail.com)
% Last Modified: November 2008
function all_neighbours = find_neighbours(data)
  % Determine X and Y coordinates of locations
  X = reshape(repmat(1:sqrt(size(data, 2)), ...
    sqrt(size(data, 2)),1),size(data, 2),1);
  Y = reshape(repmat((1:sqrt(size(data, 2)))', ...
    1,sqrt(size(data, 2))),size(data, 2),1);
  
  all_neighbours = cell(1, size(data, 2));
  
  for site = 1:size(data, 2)    
    % Find x and y coordinates of 4 direct neighbours
    x_n = zeros(1, 4); y_n = zeros(1, 4);
    
    x_n(1) = X(site); x_n(2) = X(site); 
    x_n(3) = X(site) - 1; x_n(4) = X(site) + 1; 
    
    y_n(1) = Y(site) + 1; y_n(2) = Y(site) - 1; 
    y_n(3) = Y(site); y_n(4) = Y(site); 
    
    neighbourhood = zeros(1, 4);
    
    % Find locations for those coordinates
    for i = 1:4
      if(find(X == x_n(i) & Y == y_n(i)))
        neighbourhood(i) = find(X == x_n(i) & Y == y_n(i));
      else
        neighbourhood(i) = -1;
      end
    end

    neighbourhood(neighbourhood == -1) = [];
    
    all_neighbours{site} = neighbourhood;
  end
  
  
