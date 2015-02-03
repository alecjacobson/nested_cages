function [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = matrioshka_dolls(V0,F0,levels,sep_thick,wall_thick,varargin)
  % MATRISHKA_DOLLS
  % [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = matrioshka_dolls(V0,F0,levels,sep_thick,wall_thick,varargin)
  %
  % Given a fine traingle mesh (V0,F0), an output numeber of faces
  % (for each level), separation and wall thicknesses,
  % it obtains volume minimizing Matrioshka dolls
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   levels  (#levels)x1 vector specifying the number of faces for each
  %   output mesh. We assume level(j)<level(j+1)
  %   sep_thickness   (minimum) separation between levels
  %   wall_thickness   wall thickness for each layer
  %   Optional: 
  %     'quadrature_order': 1, 2 or 3 (default=2)
  %     'beta_init': initial step size for energy minimization
  % Output:
  %   cages_V   array with (#levels) matrices with vertex positions. 
  %             cages_V{k} corresponds to levels(k)-output mesh
  %   cages_F   array with (#levels) matrices with face indices. 
  %             cages_F{k} corresponds to levels(k)-output mesh
  %   Pall      sequence of meshes from the flow
  %   V_coarse  initial coarse mesh for each level
  %   F_coarse  initial coarse mesh for each level
  %   timing:   struct with timing.decim, timing.flow and 
  %             timing.simulation
  
  
  energy_expansion = 'displacement_step';
  energy_final = 'volume';
  quadrature_order = 2;
  V_coarse = [];
  F_coarse = [];
  Pall = [];
  beta_init = 1e-2;
  step_size = 1e-3;
  % coarse mesh expansion (during flow)
  expand_every = 0;
  Fquad = [];
  partial_path = 'partial.mat';
  smoothing = 0;
  D_CV_MIN = 1e-5;
  BETA_MIN = 1e-3;

  % save timings
  timing.decimation = 0.0;
  timing.flow = 0.0;
  timing.expansion = 0.0;
  timing.etienne = 0; % not exactly timing, but let's put in this struct
  timing.final_energy = 0.0;
  debug = true;
  first_only = false;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'QuadratureOrder','StepSize','V_coarse','F_coarse','ExpandEvery', ...
      'ExpansionEnergy','FinalEnergy','BetaInit','Debug','Fquad', ...
      'PartialPath','Smoothing', ...
      'D_CV_MIN','FirstOnly','BETA_MIN'}, ...
    {'quadrature_order','step_size','V_coarse','F_coarse','expand_every', ...
      'energy_expansion','energy_final','beta_init','debug','Fquad', ...
      'partial_path','smoothing', ...
      'D_CV_MIN','first_only','BETA_MIN'});
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
  V_coarse{2*num_levels+1} = V0;
  F_coarse{2*num_levels+1} = F0;
  
  cages_V{2*num_levels+1} = V0;
  cages_F{2*num_levels+1} = F0;
  
  for k=num_levels:-1:1
      
      tic
      [V_coarse{2*k},F_coarse{2*k}] = cgal_simplification(cages_V{2*k+1},cages_F{2*k+1},levels(k));
      [V_coarse{2*k},F_coarse{2*k}] = meshfix(V_coarse{2*k},F_coarse{2*k});
      timing.decimation = timing.decimation + toc;
      
      tic
      ratio = levels(k)/size(F_coarse{2*k+1},1);
      switch char(java.lang.System.getProperty('user.name'))
      case 'ajx'
        [V_coarse{2*k},F_coarse{2*k}] = ...
          decimate_cgal(cages_V{2*k+1},cages_F{2*k+1},ratio);
      otherwise
        [V_coarse{2*k},F_coarse{2*k}] = ...
          cgal_simplification(cages_V{2*k+1},cages_F{2*k+1},levels(k));
      end
      [~,~,siIF] = selfintersect( ...
        V_coarse{2*k},F_coarse{2*k},'DetectOnly',true,'FirstOnly',true);
      if exist('meshfix','file')
        % Also remove triangles with tiny angles
        for iter = 1:100
          % 1*needed for multiplying
          A = 1*facet_adjacency_matrix(F_coarse{2*k});
          small_angles = ...
            min(internalangles(V_coarse{2*k},F_coarse{2*k}),[],2)<(5/180*pi); %5??
          if ~any(small_angles)
            break;
          end
          for grow = 1:iter-1
            small_angles = (A*small_angles)~=0;
          end
          fprintf('Removing %d skinny facets...\n',nnz(small_angles));
          F_coarse{2*k} = F_coarse{2*k}(~small_angles,:);
          [V_coarse{2*k},F_coarse{2*k}] = meshfix(V_coarse{k},F_coarse{k});
        end
      else
          error('Decimation contains self-intersections, but no meshfix');
      end
      timing.decimation = timing.decimation + toc;
      
      % flow
      tic
      [Pall,Pall_coarse] = shrink_fine_expand_coarse_3D(cages_V{2*k+1},cages_F{2*k+1},...
          V_coarse{2*k},F_coarse{2*k},'quadrature_order',quadrature_order,...
          'step_size',step_size,'expand_every',expand_every, ...
          'smoothing',smoothing,'FirstOnly',first_only);
      Pall_all_times{k} = Pall;
      timing.flow = timing.flow + toc;
      
    % push coarse mesh with physical simulation to obtain the cages
    [V_coarse_new,timing.per_layer{k},~] = ...
      combined_step_project( ...
        Pall,cages_F{2*k+1},...
        Pall_coarse(:,:,end),F_coarse{2*k}, ...
        'RefCage',V_coarse{2*k}, ...
        'ExpansionEnergy',energy_expansion, ...
        'FinalEnergy',energy_final, ...
        'EpsExpansion',sep_thick, ...
        'EpsFinal',sep_thick, ...
        'BetaInit',beta_init, ...
        'SkipElTopo',true, ...
        'Debug',debug,...
        'Fquad',Fquad,...
        'D_CV_MIN',D_CV_MIN,...
        'BETA_MIN',BETA_MIN);
    
    % output level
    cages_F{2*k} = F_coarse{2*k};
    cages_V{2*k} = V_coarse_new;
    save(partial_path,'cages_V','cages_F','timing');
    
    fprintf('done generated layer. Now generation of the wall');
    % first copy layer
    V_coarse{2*k-1} = V_coarse{2*k};
    F_coarse{2*k-1} = F_coarse{2*k};
    
    % flow
    tic
    [Pall,Pall_coarse] = shrink_fine_expand_coarse_3D(cages_V{2*k},cages_F{2*k},...
        V_coarse{2*k-1},F_coarse{2*k-1},'quadrature_order',quadrature_order,...
        'step_size',step_size,'expand_every',expand_every, ...
        'smoothing',smoothing,'FirstOnly',first_only);
    Pall_all_times{k} = Pall;
    timing.flow = timing.flow + toc;
    
    % push coarse mesh with physical simulation to obtain the cages
    [V_coarse_new,timing.per_layer{k},~] = ...
      combined_step_project( ...
        Pall,cages_F{2*k},...
        Pall_coarse(:,:,end),F_coarse{2*k-1}, ...
        'RefCage',V_coarse{2*k-1}, ...
        'ExpansionEnergy',energy_expansion, ...
        'FinalEnergy',energy_final, ...
        'EpsExpansion',wall_thick, ...
        'EpsFinal',wall_thick, ...
        'BetaInit',beta_init, ...
        'SkipElTopo',true, ...
        'Debug',debug,...
        'Fquad',Fquad,...
        'D_CV_MIN',D_CV_MIN,...
        'BETA_MIN',BETA_MIN);
    
    % output level
    cages_F{2*k-1} = F_coarse{2*k-1};
    cages_V{2*k-1} = V_coarse_new;
    save(partial_path,'cages_V','cages_F','timing');
      
  end