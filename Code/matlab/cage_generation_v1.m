function [V_coarse,F_coarse,V_shrink,F_shrink] = cage_generation_v1(V0,F0,N)
  % CAGE_GENERATION_V1
  % [V,F] = cage_generation_v1(V0,F0,N)
  %
  % Coarsens a given watertight mesh, shrinks it towards the medial
  % axis until it has no more intersections with its coarse version
  % and grows the mesh back pushing away the cage around it.
  %
  % Inputs:
  %   V0  list of surface vertex positions of exterior mesh, # vertices by 3
  %   F0  list of surface face indices of exterior triangle mesh, # faces by 3
  %   N   target number of vertices for the cage
  % Outputs:
  %   V_coarse  list of vertices of the final cage
  %   F_coarse  list of faces of the final cage
  %   V_shrink  list of vertices of the shrunk mesh
  %   F_shrink  list of faces of the shrunk mesh
  
  % coarsen the mesh with Qslim
  [~,V_coarse,F_coarse,~] = qslim(V0,F0,N);
  
  % initialize (V_shrink,F_shrink) as (V0,F0)
  V_shrink = V0;
  F_shrink = F0;
  
  % concatenate meshes and test for intersections
  V_all = [V_coarse;V_shrink];
  F_all = [F_coarse;F_shrink+size(V_coarse,1)];
  [~,~,IF] = selfintersect(V_all,F_all,'DetectOnly',true);
  % this vector  will have different values if the cage intersects
  % the shrunk mesh
  int_indices = sum(sort(IF,2)>=N);
  
  k = 1;
  while (size(IF,1)> 0 && int_indices(1)~=int_indices(2) && k<100) 
      % display message on the flow ietartion
      message = sprintf('\n \n shrinking iteration %d',k);
      disp(message);
      
      % Shrink mesh towards the medial axis
      [V_shrink,F_shrink] = mcfskel_prompt(V_shrink,F_shrink);
      k = k + 1;
      
      % concatenate meshes and test for intersections
      V_all = [V_coarse;V_shrink];
      F_all = [F_coarse;F_shrink+size(V_coarse,1)];
      [~,~,IF] = selfintersect(V_all,F_all,'DetectOnly',true);
      int_indices = sum(sort(IF,2)>=N);
      
  end
  