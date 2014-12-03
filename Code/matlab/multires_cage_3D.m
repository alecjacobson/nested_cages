function [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_cage_3D(V0,F0,levels,varargin)
  % MULTIRES_CAGE_3D
  % multires_cage_3D(V,F,varargin)
  %
  % Given a fine traingle mesh (V0,F0) and output numeber of faces )levels)
  % obtains nested cages some flow and physical simulation
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   levels  (#levels)x1 vector specifying the number of faces for each
  %   output mesh. We assume level(j)<level(j+1)
  %   Optional: 
  %     'delta_t' followed by time step for the flow and for the physical
  %     simulation
  %    'flow_type' cMCF flow used to shrink the input polygon towards
  %     inside the coarse mesh. It can be 'surface' ...
  %     'loop_collisions' followed by number of attempts to get to a
  %     collision-free state at each time step of the physical simulation
  %     'QPMethod' followed by either 'min_quad_with_fixed_active_set' or
  %       'quadprog'
  %     'energy' followed by either 'displacement_step' or
  %       'displacement_initial'
  %     'back_factor' factor that determines how much the vertices
  %     of the coarse mesh will try to go back to their initial positions
  %     when 'energy' is 'displacement_initial'
  %     ('V_coarse','F_coarse'): previously computed initial coarse layers
  %     'Pall': previously computed flow (single layer only for now)
  % Output:
  %   cages_V   array with (#levels) matrices with vertices positions. 
  %             cages_V{k} corresponds to levels(k)-output mesh
  %   cages_F   array with (#levels) matrices with face indices. 
  %             cages_F{k} corresponds to levels(k)-output mesh
  %   Pall      sequence of meshes from the flow
  %   V_coarse  initial coarse mesh
  %   F_coarse  initial coarse mesh
  
  delta_t = 0.001;
  flow_type = 'distance_steepest_descent';
  simulation_steps = 1;
  qp_method = 'quadprog';
  energy = 'displacement_step';
  loop_collisions = 50;
  back_factor = 1.0;
  simulation = 'positions';
  quadrature_order = 1;
  V_coarse = [];
  F_coarse = [];
  Pall = [];
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'delta_t'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              delta_t = varargin{ii};
          case 'flow_type'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              flow_type = varargin{ii};
          case 'simulation_steps'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              simulation_steps = varargin{ii};
          case 'energy'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              energy = varargin{ii};
          case 'loop_collisions'
            assert(ii+1<=numel(varargin));
            ii = ii+1;
            loop_collisions = varargin{ii};
          case 'back_factor'
            assert(ii+1<=numel(varargin));
            ii = ii+1;
            back_factor = varargin{ii};
          case 'QPMethod'
            ii = ii+1;
            assert(ii<=numel(varargin));
            qp_method = varargin{ii};
          case 'simulation'
            ii = ii+1;
            assert(ii<=numel(varargin));
            simulation = varargin{ii};
          case 'quadrature_order'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              quadrature_order = varargin{ii};
          case 'Pall'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              Pall = varargin{ii};
          case 'V_coarse'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              V_coarse = varargin{ii};
          case 'F_coarse'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              F_coarse = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % number of levels
  num_levels = size(levels,2);
  
  if (isempty(V_coarse)&&isempty(F_coarse))
      % coarsen the mesh with Qslim with different levels
      for k=1:num_levels
          % if energy is symmetry_x, perform special decimation to produce
          % symmetric initial coarse meshes
          if strcmp(energy,'symmetry_x')

              % generate coarse layer with qslim and fix it
              [~,V_coarse_,F_coarse_,~] = qslim(V0,F0,levels(k));
              [V_coarse_,F_coarse_] = meshfix_matlab(V_coarse_,F_coarse_);

              % intersect with x=0 plane
              [V_inter,F_inter,~] = selfintersect([V_coarse_; 0 -1 -1; 0 1 -1; 0 -1 1; 0 1 1],...
                  [F_coarse_; size(V_coarse_,1)+1 size(V_coarse_,1)+2 size(V_coarse_,1)+3;...
                  size(V_coarse_,1)+2 size(V_coarse_,1)+4 size(V_coarse_,1)+3]);

              % select left side of the mesh (x<0)
              b=1;
              while b<=size(F_inter,1)
                  if ((V_inter(F_inter(b,1),1)>=0) && (V_inter(F_inter(b,2),1)>=0) && (V_inter(F_inter(b,3),1)>=0))
                      F_inter = [F_inter(1:b-1,:);F_inter(b+1:end,:)];
                  else
                      b = b+1;
                  end
              end
              [RV,IM] = remove_unreferenced(V_inter,F_inter);
              V_inter = RV;
              F_inter = IM(F_inter);

              % duplicate and glue
              V_coarse_ = [V_inter; [-V_inter(:,1) V_inter(:,2) V_inter(:,3)]];
              F_coarse_ = [F_inter; [size(V_inter,1)+F_inter(:,1) size(V_inter,1)+F_inter(:,3) size(V_inter,1)+F_inter(:,2)]];
              [V_coarse_,~,SVJ] = remove_duplicate_vertices(V_coarse_,1e-7);
              F_coarse_ = SVJ(F_coarse_);

              V_coarse{k} = V_coarse_;
              F_coarse{k} = F_coarse_;

          else
             [~,V_coarse{k},F_coarse{k},~] = qslim(V0,F0,levels(k));
             [V_coarse{k} F_coarse{k}] = meshfix_matlab(V_coarse{k},F_coarse{k});
          end
      end
  end
  
  figure;
%   input('');
  
  % shrink meshes until they are inside the previous meshes
  % also concatenate the meshes on the way
  V_coarse_all = V_coarse{1};
  F_coarse_all = F_coarse{1};
  total_vertices = size(V_coarse_all,1);
  
  for k=2:num_levels
    %shrink
    if (isempty(Pall) || num_levels>1)
        Pall = shrink_fine_to_inside_coarse_3D(V_coarse{k},F_coarse{k},...
            V_coarse{k-1},F_coarse{k-1},'delta_t',delta_t,'flow_type',flow_type,...
            'V_to_intersect',V_coarse_all,'F_to_intersect',F_coarse_all,...
            'quadrature_order',quadrature_order);
    end
        
      % push coarse mesh with physical simulation to obtain the cages
      switch simulation
          case 'positions'
              V_coarse_new = positions_expand_fine_pushing_coarse_3d(Pall,F_coarse{k},...
              V_coarse_all,F_coarse_all,delta_t,'simulation_steps',simulation_steps,...
              'energy',energy,'loop_collisions',loop_collisions,...
              'back_factor',back_factor,'QPMethod',qp_method);
          case 'eltopo'
              [V_coarse_new,V_new] = eltopo_step_project(Pall(:,:,[1 end]),F_coarse{k},...
              V_coarse_all,F_coarse_all,delta_t,'simulation_steps',simulation_steps,'energy',energy);
          otherwise
              error('Unsupported simulation method: %s',varargin{ii});
      end
    
%       input('');
    % concatenate
    V_coarse_all = [V_new;V_coarse_new];
    V_coarse{k} = V_new;
    F_coarse_all = [F_coarse{k};size(V_coarse{k},1)+F_coarse_all];
    total_vertices = size(V_coarse_all,1);
    
    cla;
    
  end
%   meshplot(V_coarse_all,F_coarse_all);

  % shrink input mesh to inside penultimate level
  if (isempty(Pall) || num_levels>1)
      Pall = shrink_fine_to_inside_coarse_3D(V0,F0,V_coarse{num_levels},...
          F_coarse{num_levels},'delta_t',delta_t,'flow_type',flow_type,...
          'V_to_intersect',V_coarse_all,'F_to_intersect',F_coarse_all,...
          'quadrature_order',quadrature_order);
  end

  
  % push coarse mesh with physical simulation to obtain the cages
  switch simulation
      case 'positions'
          V_coarse_final = positions_expand_fine_pushing_coarse_3d(Pall,F0,...
              V_coarse_all,F_coarse_all,delta_t,'simulation_steps',simulation_steps,...
              'energy',energy,'loop_collisions',loop_collisions,...
              'back_factor',back_factor,'QPMethod',qp_method);
      case 'eltopo'
            [V_coarse_final,~] = eltopo_step_project(Pall(:,:,[1 end]),F0,...
                  V_coarse_all,F_coarse_all,delta_t,'simulation_steps',simulation_steps,'energy',energy);
          
          
  end
  
  total_vertices = size(V_coarse_final,1);
  for k=1:num_levels
      cages_V{k} = V_coarse_final(total_vertices-size(V_coarse{k},1)+1:total_vertices,:);
      cages_F{k} = F_coarse{k};
      total_vertices = total_vertices - size(V_coarse{k},1);
  end