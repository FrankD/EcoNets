% Function add_spatial_corr adds one spatial correlation node to each
% location, which corresponds to a measure of the spacial autocorrelation
% for the given species. 
%
% Autocorrelation is measured by taking the neighbours of a location, and
% taking the average of the abundance values of the given species at those 
% neighbouring locations. 
%
% Inputs:
% Data - The original data, an NxM matrix, where N is the number of species
% and M is the number of locations
% species - The number of the species for which we are adding the
% autocorrelation
% all_neighbours - A 1xM cell array containing the neighbours for each
% location. all_neighbours{i} is a matrix of size 1xL, where L is the
% number of neighbours of location i.
% 
% Outputs:
% - ext_data - The original data, extended with a new row containing the
% autocorrelation values for the given species. An N+1xM matrix.
% 
% Author: Frank Dondelinger (frank.dondelinger@gmail.com)
% Last Modified: November 2008
function ext_data = add_spatial_corr(data, species, all_neighbours)
  species_no = size(data, 1); 

  for site = 1:size(data, 2)
   neighbours = all_neighbours{site};
   
   % All neighbours are considered equally important.
   weight = 1/size(neighbours, 2);
   
   spatial_inf = 0;
   
   for neighbour = neighbours
     spatial_inf = spatial_inf + weight*data(species, neighbour);
   end
   
   data(species_no + 1, site) = spatial_inf;
  end
  
  ext_data = data;