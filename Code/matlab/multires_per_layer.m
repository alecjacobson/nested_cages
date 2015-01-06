function [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,levels,varargin)
  % MULTIRES_PER_LAYER
  % [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,levels,varargin)
  %
  % Given a fine traingle mesh (V0,F0) and an output numeber of faces
  % (for each level), obtains nested cages based on specified method
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   levels  (#levels)x1 vector specifying the number of faces for each
  %   output mesh. We assume level(j)<level(j+1)
  %   Optional: 
  %     'energy_expansion': followed by either 'displacement_step' or
  %       'displacement_initial' or 'symmetry_x' or  'volume'
  %     'energy_final': followed by either 'displacement_step' or
  %       'displacement_initial' or 'symmetry_x' or  'volume'
  %     'quadrature_order': 1, 2 or 3 (default=1)
  %     ('V_coarse','F_coarse'): previously computed initial coarse layers
  %     'eps_distance': separation attempted by the physical simulation
  %     (default=1e-4)
  %     'beta_init': step length for energy minimization (default=1e-2)
  %     'step_size': step length for the flow (default=1e-3)
  %     'expand_every': expands the coarse mesh during the flow every
  %         'expand_every' steps.
  % Output:
  %   cages_V   array with (#levels) matrices with vertices positions. 
  %             cages_V{k} corresponds to levels(k)-output mesh
  %   cages_F   array with (#levels) matrices with face indices. 
  %             cages_F{k} corresponds to levels(k)-output mesh
  %   Pall      sequence of meshes from the flow
  %   V_coarse  initial coarse mesh for each level
  %   F_coarse  initial coarse mesh for each level
  %   timing:   struct with timing.decim, timing.flow and 
  %             timing.simulation, timing.etienne
  
  energy_expansion = 'displacement_step';
  energy_final = 'volume';
  quadrature_order = 2;
  V_coarse = [];
  F_coarse = [];
  Pall = [];
  % below parameter only used for ElTopo
  eps_distance = 1e-4; 
  beta_init = 1e-2;
  step_size = 1e-3;
  % coarse mesh expansion (during flow)
  expand_every = 0;
  
  % save timings
  timing.decimation = 0.0;
  timing.flow = 0.0;
  timing.expansion = 0.0;
  timing.etienne = 0; % not exactly timing, but let's put in this struct
  timing.final_energy = 0.0;
  
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'energy_expansion'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              energy_expansion = varargin{ii};
          case 'energy_final'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              energy_final = varargin{ii};
          case 'quadrature_order'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              quadrature_order = varargin{ii};
          case 'V_coarse'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              V_coarse = varargin{ii};
          case 'F_coarse'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              F_coarse = varargin{ii};
          case 'eps_distance'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              eps_distance = varargin{ii};
          case 'beta_init'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              beta_init = varargin{ii};
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
  
    % check if layers were previously prescribed
  if (isempty(V_coarse) && isempty(F_coarse))
      decimation_given = 0;
  else
      decimation_given = 1;
  end
  
  % number of levels
  num_levels = size(levels,2);
  
  if (isempty(V_coarse)&&isempty(F_coarse))
      % coarsen the mesh with Qslim with different levels
      
      for k=1:num_levels
          % if energy is symmetry_x, perform special decimation to produce
          % symmetric initial coarse meshes
          if strcmp(energy_final,'symmetry_x')

              % I have migrated the symmetric decimation to a separate
              % function. Requires tests
              [V_coarse{k},F_coarse{k}] = symmetry_x_decimation(V0,F0,levels(k));

          end
          
      end
      
  end
  
  % Define last level as input mesh
  V_coarse{num_levels+1} = V0;
  F_coarse{num_levels+1} = F0;
  
  cages_V{num_levels+1} = V0;
  cages_F{num_levels+1} = F0;
  
  % loop over different levels
  for k=num_levels:-1:1
      
      % % only generate coarse layers if they were not prescribed
      if (~decimation_given)
          tic
          [V_coarse{k},F_coarse{k}] = cgal_simplification(cages_V{k+1},cages_F{k+1},levels(k));
          [V_coarse{k},F_coarse{k}] = meshfix(V_coarse{k},F_coarse{k});
          timing.decimation = timing.decimation + toc;
      end
      
%       % save partial result
%       save('partial.mat','Pall','V_coarse','F_coarse','V0','F0');
      
      % shirnk fine mesh, expand coarse mesh
      tic
      [Pall,Pall_coarse] = shrink_fine_expand_coarse_3D(cages_V{k+1},cages_F{k+1},...
          V_coarse{k},F_coarse{k},'quadrature_order',quadrature_order,'step_size',step_size,'expand_every',expand_every);
      timing.flow = timing.flow + toc;
      
      Pall_all_times{k} = Pall;
      
%       % save partial result
%       save('partial.mat','Pall','Pall_coarse','V_coarse','F_coarse','V0','F0');
      
      % push coarse mesh with physical simulation to obtain the cages
      tic
      [V_coarse_new,etienne_called,timing_expansion,timing_final] = combined_step_project(Pall,cages_F{k+1},...
          Pall_coarse(:,:,end),F_coarse{k},'energy_expansion',energy_expansion,'energy_final',energy_final,'eps_distance',eps_distance,'beta_init',beta_init);
      timing.expansion = timing.expansion + timing_expansion;
      
      timing.etienne = timing.etienne+etienne_called;
      timing.final_energy = timing.final_energy + timing_final;
      
      % output level
      cages_F{k} = F_coarse{k};
      cages_V{k} = V_coarse_new;
      %V_coarse{k} = cages_V{k};
      %F_coarse{k} = cages_F{k};
      
  end
  
  Pall = Pall_all_times;
          
  
  
end
