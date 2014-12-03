function cage = cage_2D(P,varargin)
  % CAGE_2D
  % cage_2D(P,varargin)
  %
  % Given a fine polygon P, obtains a cage for it using standard mean
  % curvature flow and physical simulation
  %
  % Input:
  %   P  M x 2 list of polygon vertex positions of the initial fine mesh
  %   Optional:
  %     'delta_t' followed by time step for the flow and for the physical
  %     simulation
  %     'tol_simplify'   tolerance used to build the initial coarse mesh with
  %     dpsimplify
  %     'loop_collisions' followed by number of attempts to get to a
  %     collision-free state at each time step of the physical simulation
  %     'QPMethod' followed by either 'min_quad_with_fixed_active_set' or
  %       'quadprog'
  %     'energy' followed by either 'displacement_step' or
  %       'displacement_initial'
  %     'back_factor' factor that determines how much the vertices
  %     of the coarse mesh will try to go back to their initial positions
  %     when 'energy' is 'displacement_initial'
  %     'flow_type' cMCF flow used to shrink the input polygon towards
  %     inside the coarse mesh. It can be 'curve' or 'surface' or 
  %     'medial_attraction'
  % Output:
  %   cage   Nx2 list of polygon vertex positions of the final coarse mesh
  
  delta_t = 0.25;
  tol_simplify = 0.05;
  loop_collisions = 20;
  qp_method = 'min_quad_with_fixed_active_set';
  energy = 'displacement_step';
  simulation_steps = 1;
  back_factor = 1.0;
  flow_type = 'curve';
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'delta_t'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              delta_t = varargin{ii};
          case 'tol_simplify'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              tol_simplify = varargin{ii};
          case 'loop_collisions'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              loop_collisions = varargin{ii};
          case 'QPMethod'
            ii = ii+1;
            assert(ii<=numel(varargin));
            qp_method = varargin{ii};
          case 'energy'
            ii = ii+1;
            assert(ii<=numel(varargin));
            energy = varargin{ii};
          case 'simulation_steps'
            ii = ii+1;
            assert(ii<=numel(varargin));
            simulation_steps = varargin{ii};
          case 'back_factor'
            ii = ii+1;
            assert(ii<=numel(varargin));
            back_factor = varargin{ii};
          case 'flow_type'
            ii = ii+1;
            assert(ii<=numel(varargin));
            flow_type = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % obtain initial coarse mesh
  P_coarse = dpsimplify(P,tol_simplify)
  pcoarse = plot([P_coarse(:,1)' P_coarse(1,1)],[P_coarse(:,2)' P_coarse(1,2)],'o-r','LineWidth',3);

%     % Hand-made cage for the debugging example
%     P_coarse = [-1 -1; 1 -1; 1 1; -1 1];

  % flow the fine mesh inside the coarse one until there are no more
  % collisions
  Pall = shrink_fine_to_inside_coarse(P',P_coarse','delta_t',delta_t,'flow_type',flow_type);
%   delete(pcoarse);


%     % Hand-made flow for the debugging example
%     Pall(:,:,2) = [-0.75 0.75; 0 0];
%     Pall(:,:,1) = [-1.5 1.5; 0 0];

  % push coarse mesh with physical simulation to obtain the cage
  cage = positions_expand_fine_pushing_coarse(Pall,P_coarse',delta_t,'loop_collisions',loop_collisions,...
      'QPMethod',qp_method,'energy',energy,'simulation_steps',simulation_steps,'back_factor',back_factor);
