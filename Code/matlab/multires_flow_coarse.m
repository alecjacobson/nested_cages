function [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_flow_coarse(V0,F0,levels,varargin)
  % MULTIRES_FLOW_COARSE
  % [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_flow_coarse(V0,F0,levels,varargin)
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
    {'QuadratureOrder','StepSize','ExpandEvery', ...
      'ExpansionEnergy','FinalEnergy','EpsExpansion','EpsFinal','BetaInit','Debug','Fquad', ...
      'SkipElTopo','PositiveProjection','PartialPath','Smoothing', ...
      'D_CV_MIN','FirstOnly','BETA_MIN','StepBackEvery'}, ...
    {'quadrature_order','step_size','expand_every', ...
      'energy_expansion','energy_final','eps_distance_expansion','eps_distance_final','beta_init','debug','Fquad', ...
      'skip_el_topo','positive_projection','partial_path','smoothing', ...
      'D_CV_MIN','first_only','BETA_MIN','step_back_every'});
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

  % number of levels
  num_levels = size(levels,2);

  % Define last level as input mesh
  V_coarse{num_levels+1} = V0;
  F_coarse{num_levels+1} = F0;

  cages_V{num_levels+1} = V0;
  cages_F{num_levels+1} = F0;
  timing.per_layer = cell(num_levels,1);
  % loop over different levels
  
  for k=num_levels:-1:1
      
    tic
    [~,V_coarse{k},F_coarse{k},V2W] = qslim(cages_V{k+1},cages_F{k+1},levels(k));
    timing.decimation = timing.decimation + toc;

    tic
    
    [Pall_coarse_mesh,Pall_coarse] = shrink_fine_expand_coarse_3D(V_coarse{k},F_coarse{k},...
        V_coarse{k},F_coarse{k},'quadrature_order',quadrature_order,...
        'step_size',step_size,'expand_every',expand_every, ...
        'smoothing',smoothing,'FirstOnly',first_only);
    
    Pall(:,:,1) = cages_V{k+1};
    % we are going to take 20 linear steps
    Pall(:,:,20) = Pall_coarse_mesh(V2W,:,end);
    for st=2:19
        Pall(:,:,st) = ((st-1)/19)*Pall(:,:,20)+(1-(st-1)/19)*Pall(:,:,1);
    end
        
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
        'BETA_MIN',BETA_MIN);


    % output level
    cages_F{k} = F_coarse{k};
    cages_V{k} = V_coarse_new;
    save(partial_path,'cages_V','cages_F','timing');
    
    Pall = [];
  end

end