function [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels)
  % PROGRESSIVE_HULLS_INTERSECTIONS reads all files that start with 
  % basename, computes self-intersections, interesections with the input
  % mesh and intersections between progressive hulls layers.
  %
  % [IF_self IF_input IF_levels] = progressive_hulls_intersections(basename)
  %
  % Input:
  %   basename: 'basename_0.obj' is the input layer 'basename_1.obj', ...,
  %   'basename_num_levels.obj' are the progressive hulls
  % Output:
  %   IF_input: num_levels by 1 vector such that IF_input(k) is the number 
  %      of intersections between 'basename_k.obj' and 'basename_0.obj'
  %   IF_levels: num_levels by num_levels matrix such that IF_levels(i,j)
  %      is the number of intersections between 'basename_i.obj' and 
  %      'basename_j.obj'
  
  [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
  for k=1:num_levels
    [hulls_V{k},hulls_F{k}] = load_mesh(sprintf('%s_%d.obj',basename,k));
  end
  
  IF_levels = zeros(num_levels,num_levels);
  for i=1:num_levels
      i
      for j=i:num_levels
          j
          
          V_all = [hulls_V{i};hulls_V{j}];
          F_all = [hulls_F{i};size(hulls_V{i},1)+hulls_F{j}];
          
          [SV,~,SVJ] = remove_duplicate_vertices(V_all,0.0);
          SF = SVJ(F_all);
          SF_unique = unique(sort(SF,2),'rows','stable');
          [~,~,IF] = selfintersect(SV,SF_unique);
          
          % exclude pairs of triangles that have two coplanar edges
          % (they may form a degenerate self-intersection)
          if (size(IF,1)>0)
              vol_edge_edge = [volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),2)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),3)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),2),2) SF_unique(IF(1:size(IF,1),2),3)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),2)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),3)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),2) SF_unique(IF(1:size(IF,1),2),3)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),2)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),3)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),2) SF_unique(IF(1:size(IF,1),2),3)])];
          else
              vol_edge_edge = [];
          end
          
          % exclude pairs of intersecting triangles 
          IF = IF(setdiff(1:size(IF,1),find(sum(abs(vol_edge_edge)<1e-16,2)>0)),:);
          
          
          if (size(IF,1)>0)
              vol_tets = [volume(SV,[SF_unique(IF(1:size(IF,1),1),:) SF_unique(IF(1:size(IF,1),2),1)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),:) SF_unique(IF(1:size(IF,1),2),2)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),1),:) SF_unique(IF(1:size(IF,1),2),3)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),2),:) SF_unique(IF(1:size(IF,1),1),1)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),2),:) SF_unique(IF(1:size(IF,1),1),2)])...
                  volume(SV,[SF_unique(IF(1:size(IF,1),2),:) SF_unique(IF(1:size(IF,1),1),3)])];
          else
              vol_tets = [];
          end
          
          % find intersections such that no vertex of a face lies on the
          % other face (this intersections aren't degenrate for sure)
          IF_levels(i,j) = size(find(sum(abs(vol_tets)>1e-16,2)>5),1);
          
          % Obs.: the computation below does not guarantee does not
          % necessarily find a non degenerate intersection (so I'm
          % commenting it). When a vertex of the second tirangle lie on the
          % fisrt triangle, and the oter 2 vertices define tetrahedra with
          % different signs, the second triangle may intersect the first
          % being tangent to it. So we better keep intersections that are
          % guaranteed to happen (above criterium)
%           % for the intersections such that only one vertex of a face
%           % lie on the other face chech signs of the volumes of the 
%           % other two tetrahedra
%           IF_2 = IF(find(sum(abs(vol_tets)>1e-16,2)==2),:);
%           for k=1:size(IF_2,1)
%               
%               % If first vertex of the second triangle lies on the first
%               % triangle, analyze sign of the two tetrahedra defined by the
%               % fisrt triangle and vertices 2 and 3 of the second triangle
%               if abs(volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),1)]) < 1e-16)
%                   if ((volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),2)])*volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),3)])) < 0)
%                       IF_levels(i,j) = IF_levels(i,j)+1;
%                   end
%               % If second vertex of the second triangle lies on the first
%               % triangle, analyze sign of the two tetrahedra defined by the
%               % fisrt triangle and vertices 1 and 3 of the second triangle
%               elseif abs(volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),2)]) < 1e-16)
%                   if ((volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),1)])*volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),3)])) < 0)
%                       IF_levels(i,j) = IF_levels(i,j)+1;
%                   end
%               % If third vertex of the second triangle lies on the first
%               % triangle, analyze sign of the two tetrahedra defined by the
%               % fisrt triangle and vertices 1 and 2 of the second triangle
%               elseif abs(volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),3)]) < 1e-16)
%                   if ((volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),1)])*volume(SV,[SF_unique(IF_2(k,1),:) SF_unique(IF_2(k,2),2)])) < 0)
%                       IF_levels(i,j) = IF_levels(i,j)+1;
%                   end
%               end
%               
%               % next: check if above code is correct and add "if the other"
%               % 2 vertices
%               
%           end
          
      end
  end
  % just replicate to have a symmetric matrix
  IF_levels_diag = diag(IF_levels);
  IF_levels = IF_levels + IF_levels';
  % the calculation above doubled the vaules in the diagonal,
  % replace by the original values
  for i=1:num_levels
    IF_levels(i,i) = IF_levels_diag(i);
  end
  
  IF_input = zeros(num_levels,1);
  for i=1:num_levels
      i
      V_all = [hulls_V{i};V0];
      F_all = [hulls_F{i};size(hulls_V{i},1)+F0];
      
      [SV,~,SVJ] = remove_duplicate_vertices(V_all,0.0);
      SF = SVJ(F_all);
      SF_unique = unique(sort(SF,2),'rows','stable');
      [~,~,IF] = selfintersect(SV,SF_unique);
      
      % exclude pairs of triangles that have two coplanar edges
      % (they may form a degenerate self-intersection)
      if (size(IF,1)>0)
          vol_edge_edge = [volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),2)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),3)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),2),2) SF_unique(IF(1:size(IF,1),2),3)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),2)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),3)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),1) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),2) SF_unique(IF(1:size(IF,1),2),3)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),2)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),1) SF_unique(IF(1:size(IF,1),2),3)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),2) SF_unique(IF(1:size(IF,1),1),3) SF_unique(IF(1:size(IF,1),2),2) SF_unique(IF(1:size(IF,1),2),3)])];
      else
          vol_edge_edge = [];
      end
      
      IF = IF(setdiff(1:size(IF,1),find(sum(abs(vol_edge_edge)<1e-16,2)>0)),:);
      
      
      if (size(IF,1)>0)
          vol_tets = [volume(SV,[SF_unique(IF(1:size(IF,1),1),:) SF_unique(IF(1:size(IF,1),2),1)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),:) SF_unique(IF(1:size(IF,1),2),2)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),1),:) SF_unique(IF(1:size(IF,1),2),3)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),2),:) SF_unique(IF(1:size(IF,1),1),1)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),2),:) SF_unique(IF(1:size(IF,1),1),2)])...
              volume(SV,[SF_unique(IF(1:size(IF,1),2),:) SF_unique(IF(1:size(IF,1),1),3)])];
      else
          vol_tets = [];
      end
      
      % find intersections such that no vertex of a face lies on the
      % other face (this intersections aren't degenrate for sure)
      IF_input(i) = size(find(sum(abs(vol_tets)>1e-16,2)>5),1);
      
  end