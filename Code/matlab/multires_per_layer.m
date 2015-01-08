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
  %     'QuadratureOrder': 1, 2 or 3 (default=1)
  %     ('V_coarse','F_coarse'): previously computed initial coarse layers
  %     'eps_distance': separation attempted by the physical simulation
  %     (default=1e-4)
  %     'step_size': step length for the flow (default=1e-3)
  %     'expand_every': expands the coarse mesh during the flow every
  %         'expand_every' steps.
  %     ... and any optional arguements to combined_step_project
  % Output:
  %   cages_V   array with (#levels) matrices with vertices positions. 
  %             cages_V{k} corresponds to levels(k)-output mesh
  %   cages_F   array with (#levels) matrices with face indices. 
  %             cages_F{k} corresponds to levels(k)-output mesh
  %   Pall      sequence of meshes from the flow
  %   V_coarse  initial coarse mesh for each level
  %   F_coarse  initial coarse mesh for each level
  %   timing:   struct with timing.decim, timing.flow and
  %     .per_layer  #layers list of timing structs for each layer

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
  debug = true;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'QuadratureOrder','StepSize','V_coarse','F_coarse','ExpandEvery', ...
      'ExpansionEnergy','FinalEnergy','Eps','BetaInit','Debug'}, ...
    {'quadrature_order','step_size','V_coarse','F_coarse','expand_every', ...
      'energy_expansion','energy_final','eps_distance','beta_init','debug'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
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
  timing.per_layer = cell(num_levels,1);
  % loop over different levels
  for k=num_levels:-1:1

      % % only generate coarse layers if they were not prescribed
      if (~decimation_given)
          tic
          [V_coarse{k},F_coarse{k}] = cgal_simplification(cages_V{k+1},cages_F{k+1},levels(k));
          [~,~,siIF] = selfintersect( ...
            V_coarse{k},F_coarse{k},'DetectOnly',true,'FirstOnly',true);
           if exist('meshfix','file')
            [V_coarse{k},F_coarse{k}] = meshfix(V_coarse{k},F_coarse{k});
           else
             error('Decimation contains self-intersections, but no meshfix');
           end
          timing.decimation = timing.decimation + toc;
      end

%       % save partial result
%       save('partial.mat','Pall','V_coarse','F_coarse','V0','F0');

      % shrink fine mesh, expand coarse mesh
      tic
      [Pall,Pall_coarse] = shrink_fine_expand_coarse_3D(cages_V{k+1},cages_F{k+1},...
          V_coarse{k},F_coarse{k},'quadrature_order',quadrature_order,'step_size',step_size,'expand_every',expand_every);
      timing.flow = timing.flow + toc;

      Pall_all_times{k} = Pall;

%       % save partial result
%       save('partial.mat','Pall','Pall_coarse','V_coarse','F_coarse','V0','F0');

      % push coarse mesh with physical simulation to obtain the cages
      [V_coarse_new,timing.per_layer{k}] = ...
        combined_step_project( ...
          Pall,cages_F{k+1},...
          Pall_coarse(:,:,end),F_coarse{k}, ...
          'ExpansionEnergy',energy_expansion, ...
          'FinalEnergy',energy_final, ...
          'Eps',eps_distance, ...
          'BetaInit',beta_init, ...
          'Debug',debug);


      % output level
      cages_F{k} = F_coarse{k};
      cages_V{k} = V_coarse_new;
      %V_coarse{k} = cages_V{k};
      %F_coarse{k} = cages_F{k};

  end

  Pall = Pall_all_times;
end
