function [cage_V,cage_F,Pall] = cage_3D(V0,F0,varargin)
  % CAGE_3D
  % cage_3D(V,F,varargin)
  %
  % Given a fine traingle mesh (V0,F0), obtains a cage for it using standard 
  % some flow and physical simulation
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   Optional:
  %    'faces_coarse' number of faces to be generated for the final coarse
  %    mesh (input parameter to qslim)  
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
  %     
  % Output:
  %   cage_V   (#vertices_cage)x3 list of mesh vertex positions of the 
  %   final coarse mesh
  %   cage_F   (#faces_cage)x3 list of vertex indices that form each face
  %   of the final coarse mesh
  
  faces_coarse = 100;
  delta_t = 0.001;
  flow_type = 'surface';
  simulation_steps = 1;
  qp_method = 'quadprog';
  simulation = 'positions';
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'faces_coarse'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              faces_coarse = varargin{ii};
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
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % coarsen the mesh with Qslim and fix with MeshFix
  [~,V_coarse,F_coarse,~] = qslim(V0,F0,faces_coarse);
  [V_coarse F_coarse] = meshfix_matlab(V_coarse,F_coarse);

%     % Problematic example 1
%     [V_coarse,F_coarse] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/examples/problem_Vcoarse.off');
%     % Problematic example 2
%     [V_coarse,F_coarse] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/examples/2014.05.29_problems_multires/mesh_2/problem_Vcoarse.off');
  
%   % hand-made example to debug
%   % 1) with 51 steps, this one goes through
%   V_coarse = [0.5 1 0; 1.5 1 0; 1.5 1 1];
%   F_coarse = [1 2 3];
  
  % shrink fine mesh to inside the coarse one
  Pall = shrink_fine_to_inside_coarse_3D(V0,F0,V_coarse,F_coarse,'delta_t',delta_t,'flow_type',flow_type);

%     % Problematic example 1
%     [U,G] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/examples/problem_V0.off');
%     Pall(:,:,2) = U;
%     [U,G] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/examples/problem_V1.off');
%     Pall(:,:,1) = U;
%     F0 = G;
%     % Problematic example 2
%     [U,G] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/examples/2014.05.29_problems_multires/mesh_2/problem_V0.off');
%     Pall(:,:,2) = U;
%     [U,G] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/examples/2014.05.29_problems_multires/mesh_2/problem_V1.off');
%     Pall(:,:,1) = U;
%     F0 = G;

%   % Hand-made example to debug
%   % 1) with 51 steps, this one goes through
%   F0 = [1 2 3];
%   Pall(:,:,1) = [1 0 0.1; 1 2 0.1; 1 0 1];
%   Pall(:,:,2) = [1 0 0.1; 1 0.5 0.1; 1 0 1];
  
  % push coarse mesh with physical simulation to obtain the cage
  switch simulation
      case 'positions'
          V_coarse_final = positions_expand_fine_pushing_coarse_3d(Pall,F0,V_coarse,F_coarse,delta_t,...
              'simulation_steps',simulation_steps,'energy',energy,'loop_collisions',loop_collisions,'back_factor',back_factor,'QPMethod',qp_method);
      case 'eltopo'
          V_coarse_final = eltopo_expand_fine_pushing_coarse_3d(Pall,F0,V_coarse,F_coarse,delta_t,'simulation_steps',simulation_steps);
      otherwise
          error('Unsupported simulation method: %s',varargin{ii});
  end
          
  
  % plot both meshes
  %   V_all = [V_coarse;Pall(:,:,end)];
  %   F_all = [F_coarse;F0+size(V_coarse,1)];
  % meshplot(V_all,F_all);
  
  
  cage_V = V_coarse_final;
  cage_F = F_coarse;