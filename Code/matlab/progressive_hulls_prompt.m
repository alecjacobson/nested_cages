function [hulls_V,hulls_F] = progressive_hulls_prompt(path_to_mesh,levels)
  % PROGRESSIVE_HULLS_PROMPT CallS porgessive_hulls from command line 
  % to construct a set of progressive hulls with number of faces specified
  % by the vector levels. Implemetation of the progressive hulls described
  % in Sander et al. 2001 "Silhouette Clipping"
  %
  % [hulls_V,hulls_F] = progressive_hulls_prompt(V0,F0,levels)
  %
  % Inputs:
  %   path_to_mesh  path to input mesh
  %   levels  #levels by 1 vector specifying number of faces for each layer
  % Outputs:
  %   hulls_V  array of matrices where hulls_V{1} is the vertex positions
  %   of the coarsest layer ...
  %   hulls_F  array of matrices where hulls_F{1} are vertex indices
  %   of the coarsest layer ...
  %
  % TODO: mex version
  
  path_to_progressive_hulls = '../progressive-hulls-novis/progressive_hulls_novis';
  [~,F0] = load_mesh(path_to_mesh);
  num_levels= size(levels,2);
  
  for k=size(levels,2):-1:1
      
      to_collapse = sprintf('%d',floor((size(F0,1)-levels(k))/2));
      command = [path_to_progressive_hulls ' ' path_to_mesh ' ' to_collapse];
      disp(command);
      
      [status, result] = system(command)
      if status~=0
          error(result)
      end
      
      [V_prog,F_prog] = load_mesh('progressive.obj');
      [RV,IM] = remove_unreferenced(V_prog,F_prog);
      RF = IM(F_prog);
      RF_small = RF(all(diff(sort(RF'))),:);
      
      hulls_V{num_levels+1-k} = RV;
      hulls_F{num_levels+1-k} = RF_small;
      
  end