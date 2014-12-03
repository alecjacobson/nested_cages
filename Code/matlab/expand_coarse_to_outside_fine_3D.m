function Pall = expand_coarse_to_outside_fine_3D(V0,F0,V_coarse,F_coarse,varargin)
  % EXPAND_COARSE_TO_OUTSIDE_FINE_3D
  % Pall = expand_coarse_to_outside_fine_3D(V0,F0,V_coarse,F_coarse,varargin)
  %
  % Given a COARSE triangle mesh (V_coarse,F_coarse) and a fine mesh 
  % (V0,F0), it expands (V_coarse,F_coarse) until it no longer 
  % intersects (V0,F0)
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
  %     'quadrature_order': 1, 2 or 3 (default=3)
  % Output:
  %   V_all   (#vertices)x3xsteps lists of mesh vertex positions of 
  %   the fine mesh at each time step
  
  
  quadrature_order = 3;
  % Parsing arguments
  ii = 1;
  V_to_intersect = V_coarse;
  F_to_intersect = F_coarse;
  while ii < numel(varargin)
      switch varargin{ii}
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
  V_exp = V0;
  F_exp = F0;
  
  % concatenate meshes and test for intersections
  V_all = [V_coarse;V_exp];
  F_all = [F_coarse;F_exp+size(V_coarse,1)];
  [~,~,IF] = selfintersect(V_all,F_all,'DetectOnly',true);
  
  %initialize flow
  ii = 1;
  Pall(:,:,1) = V_exp;
  
  
  % setting options for 'signed_distance_direction_flow' (the only flow
  % available for this method)
  area_initial = doublearea(V0,F0)/2;
  A_qv = grad_quadrature_to_vertices(V0,F0,area_initial,quadrature_order);
  s = 0.005;
  % needed to run 'signed_distance_direction_quadrature_matrix'
  % remove in the future
  M = massmatrix(V0,F0,'barycentric');
  w_lap = 0.0;
  L = cotmatrix(V0,F0);
      
  
  % plot general options
  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])
  
  % initialize plot handles
  hold on;
  pc = trisurf(F_to_intersect,V_to_intersect(:,1),V_to_intersect(:,2),V_to_intersect(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.05);
  pv = trisurf(F_exp,V_exp(:,1),V_exp(:,2),V_exp(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.4);
  % draw everything
  drawnow;
  
  % define plotting struct
  plot_struct = struct('V_shrink',V_exp,'F_shrink',F_exp,'V_to_intersect',V_to_intersect,...
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
      
      % the only flow we have implemented for this method is 'signed_distance_direction'
      V_exp_prev = V_exp;
      [V_exp,moving_vertices,grad_energy,quad,grad_quadrature,C,s] = signed_distance_direction_quadrature_matrix(V_exp,F_exp,1,quadrature_order,V_coarse,F_coarse,moving_vertices,'scalar_search',area_initial,A_qv,M,w_lap,L,-1);
      
      % Eltopo to remove intersections on the coarse mesh
      [V_exp,~] = collide_eltopo_mex(V_exp_prev,F_exp,V_exp,0,1e-4,1e-10);
      
      IF = intersect_other(V_to_intersect,F_to_intersect,V_exp,F_exp);
      
      plot_struct.quad_points = quad;
      plot_struct.grad_vertices = grad_energy;
      plot_struct.grad_quad = grad_quadrature;
      plot_struct.quad_closest = C;
      plot_struct.V_shrink = V_exp; 
      plot_struct.F_shrink = F_exp; 
      plot_struct.V_shrink_prev = V_exp_prev; 
      plot_struct = flow_plot_control(plot_struct,true);

%       hold off;

      
  end