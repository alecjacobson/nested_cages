function [V_coarse_final,etienne_called,time_expansion,time_final_energy]  ...
  = combined_step_project(P_all,F0,V_coarse,F_coarse,varargin)
  % COMBINED_STEP_PROJECT_3D
  % [V_coarse_final,etienne_called,time_expansion,time_final_energy]  ...
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
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
  %   coarse mesh
  %   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
  %   of the coarse mesh
  %   Optional:
  %     'simulation_steps' number of physical simulation steps to reach
  %     one step back of the flow
  %     'energy': followed by either 'displacement_step' or
  %       'displacement_initial' or 'symmetry_x' or 
  %       'displacement_initial_and_volume' (default)
  %     'min_progress': how much the mesh should evolve to continue the
  %     simulation
  % Output:
  %   V_coarse_final (#vertices_cage)x3 list of mesh vertex positions of the
  %   resulting cage
  %   etienne_called  number of times Etinne's functions are called
  %   time_final_energy  time spent in sec for final energy minimization 

  debug = true;

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
    G(sym_pairs(:,2),1) = grad_sym(sym_pairs(:,2),1) + 2*(V(sym_pairs(:,1),1)+V(sym_pairs(:,2),1));
    G(sym_pairs(:,2),2) = grad_sym(sym_pairs(:,2),2) - 2*(V(sym_pairs(:,1),2)-V(sym_pairs(:,2),2));
    G(sym_pairs(:,2),3) = grad_sym(sym_pairs(:,2),3) - 2*(V(sym_pairs(:,1),3)-V(sym_pairs(:,2),3));
    cb_data = [];
  end
  function [E,cb_data] = symmetry_x_energy(V,sym_pairs)
    G = symmetry_x_gradient(V,sym_pairs)
    E = trace(G'*G);
    cb_data = [];
  end

  function [G,cb_data] = surface_arap_gradient(V0,F,V)
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
    cb_data = [];
    [G,cb_data.E,cb_data.R] = arap_gradient(V0,F,V);
  end
  function [E,cb_data] = surface_arap_energy(V0,F,V)
    E = [];
    cb_data = [];
    [cb_data.G,E,cb_data.R] = arap_gradient(V0,F,V);
  end

  function [G,cb_data] = volumetric_arap_gradient(TV0,TT,TV)
    % Gradient of volumetric arap
    % 
    % Inputs:
    %   TV0  #TV>#V by 3 list of tet mesh vertices in initial positions: surface
    %     vertices come first TV0 = [V0;steiner]
    %   TT  #TT by 4 list of tet mesh indices into TV
    %   TV  #TV by 3  new tet mesh vertices (with surface coming first
    % Outputs:
    %   G  #V by 3 list of gradient vectors
    %   cb_data  unused output callback data field (needed to match expected
    %     prototype)

    % Compute one small step of arap flow
    % TV_flow = arap(...)
    % Double check sign and magnitude 
    cb_data = [];
  end
  function E = volumetric_arap_energy(TV0,TT,TV)
    % compute current energy of (TV,TT) w.r.t. (TV0,TT)
    cb_data = [];
  end

  function [CV_filtered, etienne_called,pc,pv] = one_step_project(V_prev,V,F,CV,CF,energy,...
          beta,eps_proximity,tol_dt,etienne_called,t,CV_orig,sym_pairs,pc,pv)
    % Given (V_prev,F) previous fine mesh, (V,F) new fine mesh and
    % (CV,CF) current coarse mesh, step from (V_prev,F) to (V,F) and
    % obtian energy minimizing (CV,CF).
    % 
    % Inputs:
    %   (V_prev,F)  previous fine mesh
    %   (V,F) next fine mesh
    %   (CV,CF) coarse mesh in the begining of the time step
    %   energy: 'displacement_step','displacement_step_and_volume', 
    %     'displacement_initial','displacement_initial_and_volume',
    %     'volume', 'symmetry_x', 'volumetric_arap'
    %   beta: initial step size for the gradient method
    %   eps_proximity: separation used for physical simulation
    %   tol_dt: tolerance for Eltopo stop trying to cut the time step
    %   etienne_called: how many times Etienne's code was called
    %   t: current flow step (for plotting purposes only)
    %   CV_orig: coarse mesh emebdding we minimize distance to 
    %   (only used for 'displacement_initial')
    %   sym_pairs: set of pairs that are symmetric w.r.t. x-axis
    %   pc: plot handle for coarse mesh (only used in debug mode)
    %   pv: plot handle for fine mesh (only used in debug mode)
    % Output:
    %   CV_filtered: new embedding for the coarse mesh
    %   etienne_called: updated Etienne
    %   pc: plot handle for coarse mesh (only used in debug mode)
    %   pv: plot handle for fine mesh (only used in debug mode)
    
     % initialize CV_filtered
     CV_filtered = CV;
    
     % for most energies call back data is not used
     cb_data = [];
     switch energy
      case {'displacement_step','displacement_step_and_volume'}
        energy_gradient = @(CV_prev,cb_data) displacement_gradient(CV_prev,CV_filtered);
        energy_value    = @(CV_prev,cb_data)   displacement_energy(CV_prev,CV_filtered);
      case {'displacement_initial','displacement_initial_and_volume'}
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
        energy_gradient = @(CV_prev,cb_data) surface_arap_gradient(CV_orig,CF,CV_prev);
        energy_value  = @(CV_prev,cb_data)     surface_arap_energy(CV_orig,CF,CV_prev);
      case {'volumetric_arap'}
        % Optimize for internal steiner points according to volumetric arap,
        % fixing surface to CV
        %
        % Note: this is done only once at beginning of this flow step. It is
        % **not** done for each call to `energy_gradient`
        %
        TV = arap(TV0,TT,1:size(CV,1),CV_filtered);
        TV_steiner = TV(size(CV,1)+1:end,:);
        energy_gradient = @(CV_prev,cb_data) volumetric_arap_gradient(TV0,TT,[CV_prev;TV_steiner],delta_t);
        energy_value  = @(CV_prev,cb_data)     volumetric_arap_energy(TV0,TT,[CV_prev;TV_steiner]);
      end

      % V_all_prev  where the meshes were
      % [V;CV]   where the meshes want to go
      % V_eltopo  where black box moved the meshes 
      % Then update CV

      % #Parameters
      %   beta  Magnitude of coarse mesh energy gradient we're attempting
      %   eps_proximity  desired separation distance (10*eps for el topo, eps
      %      for etienne)
      %   out_proximity  separation distance determined by etienne's inflation
      %

      
      % faces for all vertices of both meshes (for collision detection)
      F_all = [F;CF+size(V,1)];
      
      % initialize energy as inf
      E_val = inf;
      % Stepping in energy gradient direction until converged
      bb_iter = 1;
      beta_orig = beta;
      BETA_MIN = 1e-3;
      D_CV_MIN = 1e-5;
      %BETA_MIN = 1e-6;
      %D_CV_MIN = 1e-6;
      CV_prev = CV_filtered;
      while true
        % Update gradient on coarse mesh
        [CV_grad,cb_data] = energy_gradient(CV_prev,cb_data);

        assert(isempty(intersect_other(V_prev,F,CV_prev,CF,'FirstOnly',true)));
        [~,~,siIF] = selfintersect(CV_prev,CF,'DetectOnly',true,'FirstOnly',true);
        assert(isempty(siIF));

        % Try to take one step using el topo
%         assert(isempty(intersect_other(V_prev,F,CV_prev,CF,'FirstOnly',true)));
        [V_all_bb,rest_dt] = collide_eltopo_mex( ...
          [V_prev;CV_prev             ],F_all, ...
          [     V;CV_prev-beta*CV_grad], ...
          size(V,1),eps_proximity,tol_dt);
        % Did el topo fail?
        if (rest_dt>0.0)
          disp('ElTopo could not handle it, switching to velocityfilter')
          % Try to inflate current situation to accomodate eps_proximity
          [V_all_inf,out_proximity] = inflate_mex( ...
            [V_prev;CV_prev],F_all, ...
            size(V,1),eps_proximity);
          V_all_bb = velocity_filter_mex( ...
            V_all_inf, ... % current position (after some possible inflation)
            [V;CV_prev - beta*CV_grad], ... % desired position (after beta time step)
            F_all,size(V,1),out_proximity,0.01*out_proximity);
          etienne_called = etienne_called + 1; 
        end
        % At this point, V_all_bb contains vertex positions for both meshes
        % after the "Black Box" velocity filter (el topo + etienne)

        % Extract new positions for coarse mesh
        CV_filtered = V_all_bb(size(V,1)+1:end,:);
        % (CV_filtered,CF) should not intersect (V,F)
        assert(isempty(intersect_other(V,F,CV_filtered,CF,'FirstOnly',true)));
        [~,~,siIF] = selfintersect(CV_filtered,CF,'DetectOnly',true,'FirstOnly',true);
        assert(isempty(siIF));

        % Stop if the change in positions is tiny
        d_CV = max(normrow(CV_filtered - CV_prev));
        fprintf('d_CV:%g\n',d_CV);
        if d_CV < D_CV_MIN && bb_iter > 1
          fprintf('Max change in CV (%g) less than D_CV_MIN (%g)\n', ...
            d_CV,D_CV_MIN);
          break;
        end 

        % Compute energy at filtered positions
        E_val_prev = E_val;
        [E_val,cb_data] = energy_value(CV_filtered,cb_data);
        % Is energy decreasing (and not first run)
        if E_val < E_val_prev
          if bb_iter > 1
            % try to increase beta
            beta = min(1.1*beta,beta_orig);
          end
          fprintf('  Progress!\n');
        else
          assert(bb_iter > 1);
          % otherwise decrease beta
          beta = beta*0.5;
          % and roll back 
          CV_filtered = CV_prev;
          fprintf('  No progress...\n');
        end
%         assert(isempty(intersect_other(V,F,CV_filtered,CF,'FirstOnly',true)));

        % Stop if beta is now too small
        if beta < BETA_MIN
          fprintf('Beta (%g) less than BETA_MIN (%g)\n',beta,BETA_MIN);
          break;
        end

        if debug
          hold on;
          axis equal;
          delete(pc);
          delete(pv);
          % trisurf maintains previous axes, while tsuyrf doesn't
          pv = trisurf(F,V_prev(:,1),V_prev(:,2),V_prev(:,3),...
              'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2,'EdgeAlpha',0.2);
          pc = trisurf(CF,CV_filtered(:,1),CV_filtered(:,2),CV_filtered(:,3),...
              'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.1,'EdgeAlpha',0.2);
          title(sprintf('energy: %s, t: %d',energy,t),'FontSize',20,'Interpreter','none');
          drawnow;
          hold off;
        end

        % Updated previous positions for next iteration of loop
        CV_prev = CV_filtered;
        % Assuming we have succeeded in moving the fine mesh, then V_prev
        % should stay put
        V_prev = V;

        bb_iter = bb_iter+1;
      end
      assert(bb_iter >= 1,'must take at least one step');
      beta = beta_orig;
        
  end


  etienne_called = 0;

  energy_expansion = 'displacement_step';
  energy_final = 'volume';
  % define target cage (many times initial mesh)
  eps_distance = 1e-4;
  beta_init = 1e-2;
  % parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'energy_expansion'
              ii = ii+1;
              assert(ii<=numel(varargin));
              energy_expansion = varargin{ii};
          case 'energy_final'
              ii = ii+1;
              assert(ii<=numel(varargin));
              energy_final = varargin{ii};
          case 'eps_distance'
                assert(ii+1<=numel(varargin));
                ii = ii+1;
                eps_distance = varargin{ii};
          case 'beta_init'
                assert(ii+1<=numel(varargin));
                ii = ii+1;
                beta_init = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
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
  if strcmp(energy_final,'symmetry_x')
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
  % starting to time expansion reversing flow: we'll updated CV_filtered after
  % each reverse flow time step.
  tic;
  fprintf('Re-inflation...\n');
  for t = k-1:-1:1
    fprintf('reversing flow step %d...\n',t);
    % Previous positions of fine mesh
    V_prev = P_all(:,:,t+1);
    % Desired positions of fine mesh after step
    V = P_all(:,:,t);

    [CV_filtered,etienne_called,pc,pv] = one_step_project(V_prev,V,F,...
        CV_filtered,CF,energy_expansion,beta,eps_proximity,...
        tol_dt,etienne_called,t,CV_orig,sym_pairs,pc,pv);
  end
  time_expansion = toc;

  
  % final optimization
  fprintf('Final optimization...\n');
  tic
  [CV_filtered,etienne_called] = one_step_project(V_prev,V,F,...
        CV_filtered,CF,energy_final,beta,eps_proximity,...
        tol_dt,etienne_called,1,CV_orig,sym_pairs,pc,pv);
  time_final_energy = toc;

  V_coarse_final = CV_filtered;
  
end 
