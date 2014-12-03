function [Pall,F_refined] = shrink_fine_to_inside_coarse_3D(V0,F0,V_coarse,F_coarse,varargin)
  % SHRINK_FINE_TO_INSIDE_COARSE_3D
  % shrink_fine_to_inside_coarse_3D(V0,F0,V_coarse,F_coarse,varargin)
  %
  % Given a fine triangle mesh (V0,F0) and a simplified mesh 
  % (V_coarse,F_coarse), it shrinks (V0,F0) until it no longer 
  % intersects (V_coarse,F_coarse)
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
  %   coarse mesh
  %   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
  %   of the coarse mesh
  %   Optional:
  %     'flow_type' can be either 'surface' or ...
  %     'V_to_intersect', 'F_to_intersect': mesh to perform intersection
  %     tests at every iteration (different from V_coarse,F_coarse for 
  %     multi-layer)
  %     'quadrature_order': 1, 2 or 3 (default=3)
  % Output:
  %   V_all   (#vertices)x3xsteps lists of mesh vertex positions of 
  %   the fine mesh at each time step
  %   F_refined: new connectivity (refined version of the initial F0)
  
  
  lambda = 5e-3;
  flow_type = 'surface';
  quadrature_order = 3;
  % Parsing arguments
  ii = 1;
  V_to_intersect = V_coarse;
  F_to_intersect = F_coarse;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'flow_type'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              flow_type = varargin{ii};
          case 'V_to_intersect'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              V_to_intersect = varargin{ii};
          case 'F_to_intersect'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              F_to_intersect = varargin{ii};
          case 'quadrature_order'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              quadrature_order = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % initialize (V_shrink,F_shrink) as (V0,F0)
  V_shrink = V0;
  F_shrink = F0;
  n_fine = size(V0,1);
  
  % concatenate meshes and test for intersections
  V_all = [V_coarse;V_shrink];
  F_all = [F_coarse;F_shrink+size(V_coarse,1)];
  [~,~,IF] = selfintersect(V_all,F_all,'DetectOnly',true);
  
  %initialize flow
  ii = 1;
  Pall(:,:,1) = V_shrink;
  L = cotmatrix(V0,F0);
  
  if strcmp(flow_type,'signed_distances')
      
      [origin,distances,dx] = SDFGen_matlab(V_coarse,F_coarse,0.02,10);
      min_distances = min(min(min(distances)));
      
      vertex_distances = zeros(n_fine,1);
      vertex_weights = zeros(n_fine,1);
      
  elseif strcmp(flow_type,'distance_steepest_descent')
      
      [origin,distances,dx] = SDFGen_matlab(V_coarse,F_coarse,0.01,20);
      
  elseif strcmp(flow_type,'distance_steepest_descent_regularized')
      
      [origin,distances,dx] = SDFGen_matlab(V_coarse,F_coarse,0.01,20);
      
  elseif strcmp(flow_type,'distance_steepest_descent_tets')
      
      [SD,DV,DT,DF] = distance_on_conforming_tetmesh(V_coarse,F_coarse);

      % pre-compute gradient of the distance for each tet
      G = grad3(DV,DT);
      gradient = G*SD;
            
  elseif strcmp(flow_type,'distance_steepest_descent_regularized_tets')
      
      [SD,DV,DT,DF] = distance_on_conforming_tetmesh(V_coarse,F_coarse);
      
      % pre-compute gradient of the distance for each tet
      G = grad3(DV,DT);
      gradient = G*SD;
      
  elseif strcmp(flow_type,'signed_distance_direction')
      
      area_initial = doublearea(V0,F0)/2;
      A_qv = grad_quadrature_to_vertices(V0,F0,area_initial,quadrature_order);
%       M = sparse(1:size(V0,1),1:size(V0,1),sum(A_qv,2))
      M = massmatrix(V0,F0,'barycentric');
       
      % smoothing weight
      w_lap = 0.0;
      
      s = 0.005;
      
      
  end
  
  % plot general options
  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])
  
  % initialize plot handles
  hold on;
  pc = trisurf(F_to_intersect,V_to_intersect(:,1),V_to_intersect(:,2),V_to_intersect(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.05);
  pv = trisurf(F_shrink,V_shrink(:,1),V_shrink(:,2),V_shrink(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.4);
  % draw everything
  drawnow;
  
  % define plotting struct
  plot_struct = struct('V_shrink',V_shrink,'F_shrink',F_shrink,'V_to_intersect',V_to_intersect,...
      'F_to_intersect',F_to_intersect,'V_coarse',V_coarse,'F_coarse',F_coarse, 'quad_points',[],...
      'grad_vertices',[],'grad_quad',[],'quad_closest',[],'V_shrink_prev',V0,'IF',[],'pc',pc,'pc1',[],'pc_alpha',0.05,'pv',pv,'pv_alpha',0.4,...
      'p_quiver',[],'p_quiver_quad',[],'p_close',[],'p_quad',[],'show_weights',0);
  
  % this fixes the axis for the plots
  axis manual;
  
  % activates the camera toolbar in the figure and set principal rotation
  % axis to none
  cameratoolbar;
  cameratoolbar('SetCoordSys','none');
  
%   input('');
  hold off;
  
  
  % initialize all vertices as moving vertices
  moving_vertices = ones(size(V0,1),1);
  
  
  while (size(IF,1)> 0 && ii<100000) 
      
      ii = ii + 1
      
      switch flow_type
          case 'surface'
              % Shrink mesh towards its center with cMCF
              D = massmatrix(V_shrink,F0,'barycentric');
              S = D - lambda*L;
              V_shrink = S\(D*V_shrink);
              Pall(:,:,ii) = V_shrink;
          case 'au_skeleton'
              L = cotmatrix(V_shrink,F_shrink);
              V_shrink = [W_L*L;W_H]\[zeros(size(V0));W_H*V_shrink];
              W_L = s_L*W_L;
              areas_cur = doublearea(V_shrink,F0)/2;
              for i=1:size(V0,1)
                  face = mod(find(ismember(F0, i))+1,size(F0,1))+1;
                  W_H(i,i) = W_H_0(i,i)*sqrt(areas_0(face(1))/areas_cur(face(1)));
              end
              Pall(:,:,ii) = V_shrink;
          case 'au_skeleton_conformalized'
              V_shrink = [W_L*L;W_H]\[zeros(size(V0));W_H*V_shrink];
              W_L = s_L*W_L;
              areas_cur = doublearea(V_shrink,F0)/2;
              for i=1:size(V0,1)
                  face = mod(find(ismember(F0, i))+1,size(F0,1))+1;
                  W_H(i,i) = W_H_0(i,i)*sqrt(areas_0(face(1))/areas_cur(face(1)));
              end
              Pall(:,:,ii) = V_shrink;
          case 'mcf_skeletonization'
              [V_shrink,~] = mcfskel_prompt(V_shrink,F0);
              Pall(:,:,ii) = V_shrink;
          case 'signed_distances'
              % Shrink mesh towards its center with cMCF
              
              vertex_cells = ceil((V_shrink-ones(n_fine,1)*origin')/dx);
              for j = 1:n_fine
                vertex_distances(j) = distances(vertex_cells(j,1),vertex_cells(j,2),vertex_cells(j,3));
              end
              max(vertex_distances)
                            
              % define weights as a diagonal matrix
              for j = 1:n_fine
                  if vertex_distances(j)<-0.05
                        vertex_weights(j) = 0;
                  else
                      vertex_weights(j) = 1;
                  end
              end
              min(vertex_weights)
              W = sparse(1:n_fine,1:n_fine,vertex_weights,n_fine,n_fine);
              
%               L = cotmatrix(V_shrink,F_shrink);
              D = massmatrix(V_shrink,F0,'barycentric');
              S = D - lambda*W*L;
              V_shrink = S\(D*V_shrink);
              Pall(:,:,ii) = V_shrink;
              
          case 'distance_steepest_descent'
              
              [V_shrink,vertex_distances,grad_energy] = distance_steepest_descent(V_shrink,F0,origin,distances,dx,0.01,1,0.0,L,quadrature_order,'semi-implicit',V_coarse,F_coarse);
              Pall(:,:,ii) = V_shrink;
              
          case 'distance_steepest_descent_regularized'
              
              [V_shrink,vertex_distances,grad_energy] = distance_steepest_descent(V_shrink,F0,origin,distances,dx,0.01,1,0.002,L,quadrature_order,'semi-implicit',V_coarse,F_coarse);
              Pall(:,:,ii) = V_shrink;
              
          case 'distance_steepest_descent_tets'
              
              [V_shrink,vertex_distances,moving_vertices] = distance_steepest_descent_tets(V_shrink,F0,SD,DV,DT,0.005,1,0.0,L,quadrature_order,'semi-implicit',V_coarse,F_coarse,gradient,moving_vertices);
              Pall(:,:,ii) = V_shrink;
                            
          case 'distance_steepest_descent_regularized_tets'
              
              [V_shrink,vertex_distances,moving_vertices] = distance_steepest_descent_tets(V_shrink,F0,SD,DV,DT,0.1,quadrature_order,0.0,L,1,'semi-implicit',V_coarse,F_coarse,gradient,moving_vertices);
              Pall(:,:,ii) = V_shrink;
              
          case 'signed_distance_direction'
              
              V_shrink_prev = V_shrink;
              if (s>5e-4)
                [V_shrink,moving_vertices,grad_energy,quad,grad_quadrature,C,s] = signed_distance_direction_quadrature_matrix(V_shrink,F_shrink,1,quadrature_order,V_coarse,F_coarse,moving_vertices,'scalar_search',area_initial,A_qv,M,w_lap,L,1);
              end
              s
              w_lap
              % At every 20 iterations, refine mesh
              if (s>5e-4 && mod(ii,20)~=0)
                Pall(:,:,ii) = V_shrink;
                


              else
                  

                  % subdivide only faces that intersect with the coarse mesh 
                  % have at least one edge stretched twice the size
                  IF_fine = unique(IF(:,2));
                  int_edges = [F_shrink(IF_fine,1) F_shrink(IF_fine,2);...
                      F_shrink(IF_fine,2) F_shrink(IF_fine,3);...
                      F_shrink(IF_fine,3) F_shrink(IF_fine,1)];
                  initial_edge_lengths = normrow(Pall(int_edges(:,1),:,1)-Pall(int_edges(:,2),:,1));
                  current_edge_lengths = normrow(Pall(int_edges(:,1),:,ii-1)-Pall(int_edges(:,2),:,ii-1));
                  stretched_edges = find((current_edge_lengths./initial_edge_lengths)>1.2);
                  IF_fine = unique(IF_fine(mod(stretched_edges-1,size(IF_fine,1))+1));
                  
                  if ~isempty(IF_fine)
                      V_shrink_old = V_shrink;
                      F_shrink_old = F_shrink;
                      [V_shrink,F_shrink] = upsample(V_shrink,F_shrink,'OnlySelected',IF_fine);
                      Pall_new = zeros(size(V_shrink,1),3,ii);
                      for t=1:ii-1
                        [Pall_new(:,:,t),~] = upsample(Pall(:,:,t),F_shrink_old,'OnlySelected',IF_fine);
                      end
                      Pall = Pall_new;
                  end
                  
                  
                  disp('new number of faces')
                  size(F_shrink,1)
                  disp('new number of vertices')
                  size(V_shrink,1)
                  
%                   input('')
                  
                  % move new vertices to closest points on the surface
                  V_shrink_old = V_shrink;
                  [V_shrink] = shrink_vertices_individually(V_shrink,F_shrink,V_coarse,F_coarse,size(V_shrink_prev,1)+1:size(V_shrink,1));
                  
                  area_initial = doublearea(Pall(:,:,1),F_shrink)/2;
                  A_qv = grad_quadrature_to_vertices(V_shrink,F_shrink,area_initial,quadrature_order);
                  M = massmatrix(Pall(:,:,1),F_shrink,'barycentric');
                  L = cotmatrix(Pall(:,:,1),F_shrink);
                  
                  hold on
                  pc_new = plot3(V_shrink(size(V_shrink_prev,1)+1:size(V_shrink,1),1),V_shrink(size(V_shrink_prev,1)+1:size(V_shrink,1),2),...
                      V_shrink(size(V_shrink_prev,1)+1:size(V_shrink,1),3),'g.','markersize',10);
                  pc_new_old = plot3(V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),1),V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),2),...
                      V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),3),'r.','markersize',10);
                  p_quiver_mov = quiver3(V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),1),V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),2),...
                      V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),3),...
                        V_shrink(size(V_shrink_prev,1)+1:size(V_shrink,1),1)-V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),1),...
                        V_shrink(size(V_shrink_prev,1)+1:size(V_shrink,1),2) - V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),2),...
                        V_shrink(size(V_shrink_prev,1)+1:size(V_shrink,1),3) - V_shrink_old(size(V_shrink_prev,1)+1:size(V_shrink,1),3),0);
                  
                  delete(plot_struct.pv);
                  plot_struct.pv = trisurf(F_shrink,V_shrink(:,1),...
                      V_shrink(:,2),V_shrink(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',plot_struct.pv_alpha);
                  
                  drawnow;
                  
                  input('')
                  delete(pc_new);
                  delete(pc_new_old);
                  delete(p_quiver_mov);
                                            
                  Pall(:,:,ii) = V_shrink;
                  s = 2e-3;
                  moving_vertices = ones(size(V_shrink,1),1);

                  
              end
              
      end
      
      IF = intersect_other(V_to_intersect,F_to_intersect,V_shrink,F_shrink);
      
      plot_struct.quad_points = quad;
      plot_struct.grad_vertices = grad_energy;
      plot_struct.grad_quad = grad_quadrature;
      plot_struct.quad_closest = C;
      plot_struct.V_shrink = V_shrink; 
      plot_struct.F_shrink = F_shrink; 
      plot_struct.V_shrink_prev = V_shrink_prev; 
      plot_struct = flow_plot_control(plot_struct,true);
      
  end
  
  % update initial mesh (to be the remeshed version of it)
  V0 = Pall(:,:,1);
  F0 = F_shrink;
  F_refined = F0;