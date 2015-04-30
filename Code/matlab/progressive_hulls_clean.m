function [hulls_V,hulls_F] = progressive_hulls_clean(hulls_V,hulls_F)
  % PROGRESSIVE_HULLS_CLEAN cleans the output of progressive_hulls_mex
  % (which has a lot of duplicates)
  %
  % [hulls_V,hulls_F] = progressive_hulls_clean(hulls_V,hulls_F)
  %
  % Inputs:
  %   hulls_V  array of matrices where hulls_V{1} is the vertex positions
  %   of the coarsest layer ...
  %   hulls_F  array of matrices where hulls_F{1} are vertex indices
  %   of the coarsest layer ...
  % Outputs:
  %   hulls_V  array of matrices where hulls_V{1} is the vertex positions
  %   of the coarsest layer ...
  %   hulls_F  array of matrices where hulls_F{1} are vertex indices
  %   of the coarsest layer ...
  %
  
  num_levels = size(hulls_V,1);
  
  for k=1:num_levels
      
      [RV,IM] = remove_unreferenced(hulls_V{k},hulls_F{k});
      RF = IM(hulls_F{k});
      RF_small = RF(all(diff(sort(RF'))),:);
      
      hulls_V{k} = RV;
      hulls_F{k} = RF_small;
      
  end