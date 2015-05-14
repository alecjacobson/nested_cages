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
  %     ... and any optional arguements to combined_step_project or
  %       shrink_fine_expand_coarse_3D
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
  eps_distance_final = 1e-4;
  eps_distance_expansion = eps_distance_final;
  beta_init = 1e-2;
  step_size = 1e-3;
  % coarse mesh expansion (during flow)
  expand_every = 0;
  Fquad = [];
  partial_path = 'partial.mat';
  smoothing = 0;
  D_CV_MIN = 1e-5;
  BETA_MIN = 1e-3;
  step_back_every = 1;
  pflow = inf;
  flow_step_vis = false;
  shrink_points = false;

  % save timings
  timing.decimation = 0.0;
  timing.flow = 0.0;
  timing.expansion = 0.0;
  timing.etienne = 0; % not exactly timing, but let's put in this struct
  timing.final_energy = 0.0;
  debug = true;
  skip_el_topo = false;
  positive_projection = false;
  first_only = false;
  new_flow = false;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'QuadratureOrder','StepSize','V_coarse','F_coarse','ExpandEvery', ...
      'ExpansionEnergy','FinalEnergy','EpsExpansion','EpsFinal','BetaInit','Debug','Fquad', ...
      'SkipElTopo','PositiveProjection','PartialPath','Smoothing', ...
      'D_CV_MIN','FirstOnly','NewFlow','BETA_MIN','StepBackEvery','pflow','FlowStepVis','ShrinkPoints'}, ...
    {'quadrature_order','step_size','V_coarse','F_coarse','expand_every', ...
      'energy_expansion','energy_final','eps_distance_expansion','eps_distance_final','beta_init','debug','Fquad', ...
      'skip_el_topo','positive_projection','partial_path','smoothing', ...
      'D_CV_MIN','first_only','new_flow','BETA_MIN','step_back_every','pflow','flow_step_vis','shrink_points'});
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
          if (strcmp(energy_final,'symmetry_x')||strcmp(energy_expansion,'symmetry_x'))

              % I have migrated the symmetric decimation to a separate
              % function. Requires tests
              [V_coarse{k},F_coarse{k}] = symmetry_x_decimation(V0,F0,levels(k));
              decimation_given = true;

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
      ratio = levels(k)/size(F_coarse{k+1},1);
      switch char(java.lang.System.getProperty('user.name'))
      case 'ajx'
        [V_coarse{k},F_coarse{k}] = ...
          decimate_cgal(cages_V{k+1},cages_F{k+1},ratio);
      otherwise
        [V_coarse{k},F_coarse{k}] = ...
          cgal_simplification(cages_V{k+1},cages_F{k+1},levels(k));
      end
      [~,~,siIF] = selfintersect( ...
        V_coarse{k},F_coarse{k},'DetectOnly',true,'FirstOnly',true);
      if exist('meshfix','file')
        % Also remove triangles with tiny angles
%         for iter = 1:100
%           % 1*needed for multiplying
%           A = 1*facet_adjacency_matrix(F_coarse{k});
%           small_angles = ...
%             min(internalangles(V_coarse{k},F_coarse{k}),[],2)<(5/180*pi); %5??
%           if ~any(small_angles)
%             break;
%           end
%           for grow = 1:iter-1
%             small_angles = (A*small_angles)~=0;
%           end
%           fprintf('Removing %d skinny facets...\n',nnz(small_angles));
%           F_coarse{k} = F_coarse{k}(~small_angles,:);
          [V_coarse{k},F_coarse{k}] = meshfix(V_coarse{k},F_coarse{k});
%         end
      else
        error('Decimation contains self-intersections, but no meshfix');
      end
      timing.decimation = timing.decimation + toc;
    end

    tic
%     if (~new_flow)
%         % shrink fine mesh, expand coarse mesh
        [Pall,Pall_coarse] = shrink_fine_expand_coarse_3D(cages_V{k+1},cages_F{k+1},...
            V_coarse{k},F_coarse{k},'quadrature_order',quadrature_order,...
            'step_size',step_size,'expand_every',expand_every, ...
            'smoothing',smoothing,'FirstOnly',first_only,'pflow',pflow,...
            'FlowStepVis',flow_step_vis,'ShrinkPoints',shrink_points);
%     else
%         [Pall,Pall_coarse] = shrink_fine_expand_coarse_3D(V_coarse{k+1},F_coarse{k+1},...
%             V_coarse{k},F_coarse{k},'quadrature_order',quadrature_order,...
%             'step_size',step_size,'expand_every',expand_every, ...
%             'smoothing',smoothing,'FirstOnly',first_only);
%         if (k~=num_levels)
%             % shift flow steps forward
%             Pall(:,:,size(Pall_expansion,3)+1:size(Pall_expansion,3)+size(Pall,3)) = Pall;
%             % Append steps comming from previous attempt
%             Pall(:,:,1:size(Pall_expansion,3)) = Pall_expansion; 
%         end
%     end
    frames = 1:step_back_every:(size(Pall,3)-1);
    frames = [frames size(Pall,3)];
    Pall = Pall(:,:,frames);
    Pall_all_times{k} = Pall;
    timing.flow = timing.flow + toc;

    %% save partial result
    %save('partial.mat','Pall','Pall_coarse','V_coarse','F_coarse','V0','F0');

    % push coarse mesh with physical simulation to obtain the cages
    [V_coarse_new,timing.per_layer{k},~] = ...
      combined_step_project( ...
        Pall,cages_F{k+1},...
        Pall_coarse(:,:,end),F_coarse{k}, ...
        'RefCage',V_coarse{k}, ...
        'ExpansionEnergy',energy_expansion, ...
        'FinalEnergy',energy_final, ...
        'EpsExpansion',eps_distance_expansion, ...
        'EpsFinal',eps_distance_final, ...
        'BetaInit',beta_init, ...
        'SkipElTopo',skip_el_topo, ...
        'Debug',debug,...
        'Fquad',Fquad,...
        'D_CV_MIN',D_CV_MIN,...
        'BETA_MIN',BETA_MIN,...
        'ShrinkPoints',shrink_points);


    % output level
    cages_F{k} = F_coarse{k};
    cages_V{k} = V_coarse_new;
    save(partial_path,'cages_V','cages_F','timing');
  end

  Pall = Pall_all_times;
end
