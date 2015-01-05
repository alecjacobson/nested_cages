function [Pall_fine,Pall_coarse,F_exp,F_shrink] = shrink_fine_expand_coarse_3D(V0,F0,V_coarse,F_coarse,varargin)
  % SHRINK_FINE_EXPAND_COARSE_3D
  % [Pall_fine,Pall_coarse] = shrink_fine_expand_coarse_3D(V0,F0,V_coarse,F_coarse,varargin)
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
  %     'quadrature_order': 1, 2 or 3 (default=2)
  %     'step_size': time step for the flow (default=1e-3)
  %     'eps_distance': sepration bewteen shrunk mesh and coarse mesh
  %     (ongoing)
  %     'expand_every': expand coarse mesh every eexapnd_every' flow steps
  % Output:
  %   Pall_fine   (#vertices)x3xsteps lists of mesh vertex positions of 
  %   the fine mesh at each time step
  %   Pall_coarse   (#vertices_cage)x3xsteps lists of mesh vertex positions of 
  %   the fine mesh at each time step
  
  % plot general options
  cla;
  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])
  hold on
  cameratoolbar;
  cameratoolbar('SetCoordSys','none');

  quadrature_order = 2;
  eps_distance = 1e-4;
  step_size = 1e-3;
  expand_every = 0;
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'quadrature_order'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              quadrature_order = varargin{ii};
          case 'eps_distance'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              eps_distance = varargin{ii};
          case 'step_size'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              step_size = varargin{ii};
          case 'expand_every'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              expand_every = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % initialize (V_exp,F_exp) as (V_coarse,F_coarse)
  V_exp = V_coarse;
  F_exp = F_coarse;
  
%   tsurf(F_exp,V_exp,'FaceIndices',true,'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.5);
  
  % initialize (V_shrink,F_shrink) as (V0,F0)
  V_shrink = V0;
  F_shrink = F0;
  
  % concatenate meshes and test for intersections
  V_all = [V_exp;V_shrink];
  F_all = [F_exp;F_shrink+size(V_exp,1)];
  [~,~,IF] = selfintersect(V_all,F_all,'DetectOnly',true);
  
  %initialize flow
  ii = 1;
  Pall_coarse(:,:,1) = V_exp;
  Pall_fine(:,:,1) = V_shrink;
  
  
  % flow settings for the coarse mesh
  area_initial_coarse = doublearea(V_coarse,F_coarse)/2;
  A_qv_coarse = grad_quadrature_to_vertices(V_coarse,F_coarse,area_initial_coarse,quadrature_order);
  % needed to run 'signed_distance_direction_quadrature_matrix'
  % remove in the future
  M_coarse = massmatrix(V_coarse,F_coarse,'barycentric');
  w_lap_coarse = 0.0;
  L_coarse = cotmatrix(V_coarse,F_coarse);
  % initialize all vertices as moving vertices
  moving_vertices_coarse = ones(size(V_exp,1),1);
  
  % flow settings for the fine mesh
  area_initial_shrink = doublearea(V_shrink,F_shrink)/2;
  A_qv_shrink = grad_quadrature_to_vertices(V_shrink,F_shrink,area_initial_shrink,quadrature_order);
  % needed to run 'signed_distance_direction_quadrature_matrix'
  % remove in the future
  M_shrink = massmatrix(V_shrink,F_shrink,'barycentric');
  w_lap_shrink = 0.0;
  L_shrink = cotmatrix(V_shrink,F_shrink);
  % initialize all vertices as moving vertices
  moving_vertices_shrink = ones(size(V_shrink,1),1);
  
  % plot general options
  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])
  
  % initialize plot handles
  hold on;
  pc = trisurf(F_exp,V_exp(:,1),V_exp(:,2),V_exp(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.05);
  pv = trisurf(F_shrink,V_shrink(:,1),V_shrink(:,2),V_shrink(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.4);
  % draw everything
  drawnow;
  
  % define plotting struct
  plot_struct = struct('V_shrink',V_shrink,'F_shrink',F_shrink,'V_to_intersect',V_exp,...
      'F_to_intersect',F_exp,'V_coarse',V_exp,'F_coarse',F_exp, 'quad_points',[],...
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
  % for now, let's remove distance checking (some models were not going 
  % the distance inside)
  dist_int = 1e+16;
  
  while (((size(IF,1)>0 || dist_int<eps_distance)) || winding_number(V_coarse,F_coarse,V_shrink(1,:))<1e-10)
     
          ii = ii + 1;

          % first one step to expand coarse mesh
          V_exp_prev = V_exp;

          if (mod(ii,expand_every)==0)
              [V_exp,moving_vertices_coarse,grad_energy,quad,grad_quadrature,C,~] = signed_distance_direction_quadrature_matrix(V_exp,F_exp,1,quadrature_order,...
                  V_coarse,F_coarse,moving_vertices_coarse,A_qv_coarse,M_coarse,w_lap_coarse,L_coarse,-1);

              % Eltopo to remove intersections on the coarse mesh
              [V_exp,~] = collide_eltopo_mex(V_exp_prev,F_exp,V_exp,0,1e-4,1e-10);

              Pall_coarse(:,:,ii) = V_exp;

              IF = intersect_other(V_shrink,F_shrink,V_exp,F_exp);
              fprintf('%05d: number of intersections after coarse mesh expansion: %d \n',ii,size(IF,1));

              plot_struct.quad_points = quad;
              plot_struct.grad_vertices = grad_energy;
              plot_struct.grad_quad = grad_quadrature;
              plot_struct.quad_closest = C;
              plot_struct.V_coarse = V_exp;
              plot_struct.F_coarse = F_exp;
              plot_struct.V_to_intersect = V_exp;
              plot_struct.F_to_intersect = F_exp;
              %       plot_struct.V_shrink_prev = V_exp_prev;
              plot_struct = flow_plot_control(plot_struct,false);


          else

              V_shrink_prev = V_shrink;
              [V_shrink,moving_vertices_shrink,grad_energy,quad,grad_quadrature,C,~] = signed_distance_direction_quadrature_matrix(V_shrink,F_shrink,1,quadrature_order,...
                  V_exp,F_coarse,moving_vertices_shrink,A_qv_shrink,M_shrink,w_lap_shrink,L_shrink,1,'step',step_size);

              plot_struct.quad_points = quad;
              plot_struct.grad_vertices = grad_energy;
              plot_struct.grad_quad = grad_quadrature;
              plot_struct.quad_closest = C;
              plot_struct.V_shrink = V_shrink;
              plot_struct.V_shrink_prev = V_shrink_prev;
              plot_struct.F_shrink = F_shrink;
              plot_struct = flow_plot_control(plot_struct,false);

              Pall_fine(:,:,end+1) = V_shrink;
              
              V_shrink_prev = V_shrink;
              IF = intersect_other(Pall_fine(:,:,end),F_shrink,Pall_coarse(:,:,end),F_coarse);
              fprintf('%05d: number of intersections after fine mesh shrinking: %d \n',ii,size(IF,1));
%               if (size(IF,1)==0)
%                   V_int = [V_shrink;Pall_coarse(:,:,end)];
%                   F_int = [F_shrink; size(V_shrink,1)+F_coarse];
%                   dist_int = self_distance_mex(V_int,F_int,size(V_shrink,1));
%                   fprintf('distance between coarse and fine mesh %g \n',dist_int);
%                   input('');
%               end

          end

    end
  
  end
