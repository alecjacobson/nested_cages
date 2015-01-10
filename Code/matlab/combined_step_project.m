function [V_coarse_final,timing] = ...
  combined_step_project(P_all,F0,V_coarse,F_coarse,varargin)
  % COMBINED_STEP_PROJECT_3D
  % [V_coarse_final,timing]  ...
  % = combined_step_project(P_all,F0,V_coarse,F_coarse,varargin)
  %
  % Given the evolution of a fine mesh V (the last is assumed to be
  % inside the coarse mesh), grow the fine mesh to its initial position
  % pushing (V_coarse,F_coarse) in a way that results in no collisions,
  % using the ElTopo library and if Eltopo can't handle it, use Eienne's
  % velocityfilter.
  %
  % Input:
  %   P_all  (#vertices)x3xsteps list of mesh vertex positions of the initial
  %     fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the initial
  %     mesh
  %   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the coarse
  %     mesh
  %   F_coarse   (#faces_cage)x3 list of vertex indices that form each face of
  %     the coarse mesh
  %   Optional:
  %     'ExpansionEnergy': followed by either 
  %       'displacement_step'  "Drag energy" tries to stay put at each iteration
  %       'displacement_path'  "Drag energy" tries to not move at all
  %         throughout iteration: only takes one step.
  %       'displacement_initial' Tries to return to intial position L2 norm
  %       'surface_arap'  Minimize ARAP spokes and rims energy on surface
  %       'volumetric_arap'  Minimize ARAP in tet mesh inside shape
  %     'FinalEnergy' followed by energy to use for final energy after 
  %       expansion. Same options as 'ExpansionEnergy'
  %      ... and any optional arguments to one_step_project
  % Output:
  %   V_coarse_final (#vertices_cage)x3 list of mesh vertex positions of the
  %   resulting cage
  %   timing  struct containing timing information
  %     .per_step  #steps list of timing information structs per layer
  %     .expansion  total time for expansion
  %     .final  total time for final optimization
  %

  % Computes area-weighted (unnormalized) sum of face-normals incident on each
  % vertex:
  %
  % Inputs:
  %   CV  #CV by 3 list of mesh vertex positions
  %   CF  #CF by 3 list of mesh triangle indices into CV
  % Outputs:
  %   grad_vol  #CV by 3 area-weight normal sums
  function grad_vol = area_weighted_normal(CV,CF)
    % "un-normalized" normals are in fact unit normals times twice the faces area
    % (result of cross product) thus we just need to divide by 2 here
    N = normals(CV,CF)/2;
    grad_vol = full(sparse( ...
      repmat(CF(:),1,3),repmat(1:3,numel(CF),1),repmat(N,3,1),size(CV,1),3));
  end

  function [G,cb_data] = displacement_gradient(V,V0)
    % Gradient of total displacement from (V0,F) to (V,F)
    %
    % Inputs:
    %   V  #V by 3 current positions
    %   V0  #V by 3 positions to take displacement with respect to (e.g. initial
    %     positions)
    % Outputs:
    %   G  #V by 3 list of gradient vectors
    %   cb_data  unused output callback data field (needed to match expected
    %     prototype)
    G = V0-V;
    cb_data = [];
  end
  function [E,cb_data] = displacement_energy(V,V0)
    G = displacement_gradient(V,V0);
    E = trace(G'*G);
    cb_data = [];
  end

  function [G,cb_data] = volume_gradient(V,F)
    % Gradient of volume of mesh (V,F)
    %
    % Inputs:
    %   V  #V by 3 mesh positions
    %   F  #F by 3 mesh facet indices
    % Outputs:
    %   G  #V by 3 list of gradient vectors
    %   cb_data  unused output callback data field (needed to match expected
    %     prototype)
    G = area_weighted_normal(V,F);
    cb_data = [];
  end
  function [E,cb_data] = volume_energy(V,F)
    [~,E] = centroid(V,F);
    cb_data = [];
  end
  function [G,cb_data] = symmetry_x_gradient(V,sym_pairs)
    % Gradient of energy compairing differences of pairs of vertices (after
    % reflection due to symmetry over yz-plane)
    %
    % Inputs:
    %   V  #V by 3 current positions
    %   sym_pairs  #sym_pairs by 2 list of indices into V of symmetric pairs
    % Outputs:
    %   G  #V by 3 list of gradient vectors
    %   cb_data  unused output callback data field (needed to match expected 
    %     prototype)
    G = zeros(size(V));
    G(sym_pairs(:,1),1) = 2*(V(sym_pairs(:,1),1)+V(sym_pairs(:,2),1));
    G(sym_pairs(:,1),2) = 2*(V(sym_pairs(:,1),2)-V(sym_pairs(:,2),2));
    G(sym_pairs(:,1),3) = 2*(V(sym_pairs(:,1),3)-V(sym_pairs(:,2),3));
    G(sym_pairs(:,2),1) = G(sym_pairs(:,2),1) + 2*(V(sym_pairs(:,1),1)+V(sym_pairs(:,2),1));
    G(sym_pairs(:,2),2) = G(sym_pairs(:,2),2) - 2*(V(sym_pairs(:,1),2)-V(sym_pairs(:,2),2));
    G(sym_pairs(:,2),3) = G(sym_pairs(:,2),3) - 2*(V(sym_pairs(:,1),3)-V(sym_pairs(:,2),3));
    cb_data = [];
  end
  function [E,cb_data] = symmetry_x_energy(V,sym_pairs)
    G = symmetry_x_gradient(V,sym_pairs);
    E = trace(G'*G);
    cb_data = [];
  end

  function [G,cb_data] = surface_arap_gradient(V0,F,V,arap_data)
    % Gradient of surface based arap energy
    %
    % Inputs:
    %   V0  #V by 3 list of rest-pose mesh vertex positions
    %   F  #F by 3 list of mesh facet indices
    %   V  #V by 3 list of deformed mesh vertex positions
    % Outputs
    %   G  #V by 3 list of gradient vectors
    %   cb_data  callback data set with energy data.E and best fit rotations
    %     data.R
    %
    cb_data.arap_data = arap_data;
    [G,cb_data.E,cb_data.R,cb_data.arap_data] = ...
      arap_gradient(V0,F,V,'Data',cb_data.arap_data,'SinglePrecision',true);
    %warning('Saving arap state');
    %R = cb_data.R;
    %save('arap.mat','V0','F','V','G','R');
  end
  function [E,cb_data] = surface_arap_energy(V0,F,V,arap_data)
    [G,cb_data] = surface_arap_gradient(V0,F,V,arap_data);
    E = cb_data.E;
  end

  function [G,cb_data] = surface_arap_planarity_gradient(V0,F,V,arap_data,Fquad)
    cb_data.arap_data = arap_data;
    [G_arap,cb_data.E,cb_data.R,cb_data.arap_data] = ...
      arap_gradient(V0,F,V,'Data',cb_data.arap_data);
    [G_plan,~] = planarity_gradient(V,Fquad);
    G = G_arap+G_plan;
  end
  function [E,cb_data] = surface_arap_planarity_energy(V0,F,V,arap_data,Fquad)
    [G,cb_data] = surface_arap_gradient(V0,F,V,arap_data);
    E = cb_data.E;
    [E_plan,~] = planarity_energy(V,Fquad);
    E = E+E_plan;
  end

  function [G,cb_data] = volumetric_arap_gradient(V,cb_data)
    % Gradient of volumetric arap
    % 
    % Inputs:
    %   V  #V by 3 list of deformed mesh vertex positions
    %   cb_data  struct containing
    %     .TV0  #TV>#V by 3 list of tet mesh vertices in initial positions: surface
    %       vertices come first TV0 = [V0;steiner]
    %     .TT  #TT by 4 list of tet mesh indices into TV
    %     .TV  #TV by 3  new tet mesh vertices (with surface coming first
    % Outputs:
    %   G  #V by 3 list of gradient vectors
    %   cb_data  unused output callback data field (needed to match expected
    %     prototype)
    b = 1:size(V,1);
    % Save solution as initial value for next call
    [cb_data.TV,cb_data.arap_data,~,R] = arap( ...
      cb_data.TV0,cb_data.TT,b,V,'Energy','elements','V0',cb_data.TV,'Tol',1e-7,'MaxIter',inf, ...
      'Data',cb_data.arap_data);
    [G,cb_data.E] = arap_gradient(cb_data.TV0,cb_data.TT,cb_data.TV, ...
      'Energy','elements','Rotations',R,'Data',cb_data.arap_data);
    G = G(b,:);
  end
  function [E,cb_data] = volumetric_arap_energy(V,cb_data);
    % compute current energy of (TV,TT) w.r.t. (TV0,TT)
    TV_prev = cb_data.TV;
    [~,cb_data] = volumetric_arap_gradient(V,cb_data);
    % Just pass along initial values (calls to energy never change solution)
    cb_data.TV = TV_prev;
    E = cb_data.E;
  end

  function [energy_gradient,energy_value,cb_data] = ...
    energy_handles_from_string(energy)
    % for most energies call back data is not used
    cb_data = [];
    switch energy
      case {'displacement_step'}
        % Displacement to solution from previous step
        energy_gradient = @(CV_prev,cb_data) displacement_gradient(CV_prev,CV_filtered);
        energy_value    = @(CV_prev,cb_data)   displacement_energy(CV_prev,CV_filtered);
      case {'displacement_path'}
        % Displacement to self: 0 (lazy, should just return 0s but this costs
        % nothing)
        energy_gradient = @(CV_prev,cb_data) displacement_gradient(CV_prev,CV_prev);
        energy_value    = @(CV_prev,cb_data)   displacement_energy(CV_prev,CV_prev);
      case {'displacement_initial'}
        energy_gradient = @(CV_prev,cb_data) displacement_gradient(CV_prev,CV_orig);
        energy_value    = @(CV_prev,cb_data)   displacement_energy(CV_prev,CV_orig);
      case {'volume'}
        energy_gradient = @(CV_prev,cb_data) volume_gradient(CV_prev,CF);
        energy_value    = @(CV_prev,cb_data)   volume_energy(CV_prev,CF);
      case {'symmetry_x'}
        energy_gradient = @(CV_prev,cb_data) symmetry_x_gradient(CV_prev,sym_pairs);
        energy_value    = @(CV_prev,cb_data)   symmetry_x_energy(CV_prev,sym_pairs);
      case {'surface_arap'}
        % Call back data will be used to reduce calls to `fit_rotations`: this
        cb_data.arap_data = [];
        energy_gradient = @(CV_prev,cb_data) ...
          surface_arap_gradient(CV_orig,CF,CV_prev,cb_data.arap_data);
        energy_value  = @(CV_prev,cb_data) ...
          surface_arap_energy(CV_orig,CF,CV_prev,cb_data.arap_data);
      case {'volumetric_arap'}
        % Optimize for internal steiner points according to volumetric arap,
        % fixing surface to CV
        %
        % Note: this is done only once at beginning of this flow step. It is
        % **not** done for each call to `energy_gradient`
        %
        [cb_data.TV0,cb_data.TT,cb_data.TF] = tetgen(CV_orig,CF,'Flags','-q2Y');
        cb_data.TV = cb_data.TV0;
        cb_data.arap_data = [];
        energy_gradient = @(CV_prev,cb_data) volumetric_arap_gradient(CV_prev,cb_data);
        energy_value  = @(CV_prev,cb_data) volumetric_arap_energy(CV_prev,cb_data);
      case {'planarity'}
        energy_gradient = @(CV_prev,cb_data) planarity_gradient(CV_prev,Fquad);
        energy_value    = @(CV_prev,cb_data)   planarity_energy(CV_prev,Fquad);
      case {'surface_arap_planarity'}
        cb_data.arap_data = [];
        energy_gradient = @(CV_prev,cb_data) surface_arap_planarity_gradient(CV_orig,CF,CV_prev,cb_data.arap_data,Fquad);
        energy_value    = @(CV_prev,cb_data) surface_arap_planarity_energy(CV_orig,CF,CV_prev,cb_data.arap_data,Fquad);
      end
    end

  % default values
  energy_expansion = 'displacement_step';
  energy_final = 'volume';
  % define target cage (many times initial mesh)
  eps_distance = 1e-4;
  beta_init = 1e-2;
  debug = true;
  Fquad = [];
  skip_el_topo = false;
  D_CV_MIN = 1e-5;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {  'SkipElTopo','ExpansionEnergy','FinalEnergy','Eps','BetaInit','Debug','Fquad','D_CV_MIN'}, ...
    {'skip_el_topo','energy_expansion','energy_final','eps_distance','beta_init','debug','Fquad','D_CV_MIN'});
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

  % fine mesh
  F = F0;
  % coarse  mesh
  CV = V_coarse;
  CF = F_coarse;
  % Shrunk fine mesh
  V = P_all(:,:,end);

  % activates the camera toolbar in the figure and set principal rotation
  % axis to none
  cameratoolbar;
  cameratoolbar('SetCoordSys','none');

  % initialize plot handles
  pc = [];
  pv = [];

  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])

  cla;

  % number of steps
  k = size(P_all,3);

  % initialize beta
  beta = beta_init;

  % Eltopo parameters
  eps_proximity = eps_distance;
  tol_dt = 1e-1;

  sym_pairs = [];
  if (strcmp(energy_final,'symmetry_x')||strcmp(energy_expansion,'symmetry_x'))
      % brute-force way of detecting pairs of points symmetric w.r.t.
      % to x=0 plane
      for i=1:size(CV,1)
          p = CV(i,:);
          for j=i:size(CV,1)
              if (abs(p(1)+CV(j,1))<1e-4 && abs(p(2)-CV(j,2))<1e-4 && abs(p(3)-CV(j,3))<1e-24)
                  sym_pairs = [sym_pairs; i j];
              end
          end
      end
  end

  % Initial guess at filtered coarse mesh vertices
  CV_orig = CV;
  CV_filtered = CV;
  %  Important that this is called after CV_orig and sym_pairs is set
  [energy_gradient,energy_value,cb_data] = ...
    energy_handles_from_string(energy_expansion);
  % starting to time expansion reversing flow: we'll updated CV_filtered after
  % each reverse flow time step.
  tic;
  timing.per_step = cell(k,1);
  fprintf('Re-inflation...\n');
  for t = k-1:-1:1
    fprintf('reversing flow step %d...\n',t);
    plot_info.t = t;
    plot_info.energy = energy_expansion;
    % Previous positions of fine mesh
    V_prev = P_all(:,:,t+1);
    % Desired positions of fine mesh after step
    V = P_all(:,:,t);
    if ~strcmp(energy_expansion,'displacement_path')
        [CV_filtered,timing.per_step{k-t}] = one_step_project( ...
          V_prev,V,F,...
          CV_filtered,CF, ...
          energy_gradient,energy_value,cb_data, ...
          'BetaInit',beta,'Eps',eps_proximity,...
          'SkipElTopo',skip_el_topo, ...
          'Tol',tol_dt,'PlotInfo',plot_info,'Debug',debug,'D_CV_MIN',D_CV_MIN);
    else
        disp('Re-inflating with 10*eps')
        [CV_filtered,timing.per_step{k-t}] = one_step_project( ...
          V_prev,V,F,...
          CV_filtered,CF, ...
          energy_gradient,energy_value,cb_data, ...
          'BetaInit',beta,'Eps',10*eps_proximity,...
          'SkipElTopo',skip_el_topo, ...
          'Tol',tol_dt,'PlotInfo',plot_info,'Debug',debug,'D_CV_MIN',D_CV_MIN);
    end
  end
  timing.expansion = toc;


  % final optimization
  tic
  if ~strcmp(energy_final,'none')
    plot_info.t = inf;
    plot_info.energy = energy_final;
    fprintf('Final optimization...\n');
    [energy_gradient,energy_value,cb_data] = ...
      energy_handles_from_string(energy_final);
    [CV_filtered,timing.per_step{k}] = ...
      one_step_project( ...
        V,V,F, ...
        CV_filtered,CF, ...
        energy_gradient,energy_value,cb_data, ...
        'SkipElTopo',skip_el_topo, ...
        'BetaInit',beta,'Eps',eps_proximity, ...
        'Tol',tol_dt,'PlotInfo',plot_info,'Debug',debug,'D_CV_MIN',D_CV_MIN);
  end
  timing.final_energy = toc;

  V_coarse_final = CV_filtered;

end 
