function [V_coarse_final,V_new,F_to_refine,etienne_called,time_expansion,time_final_energy]  ...
  = combined_step_project(P_all,F0,V_coarse,F_coarse,varargin)
  % COMBINED_STEP_PROJECT_3D
  % [V_coarse_final,V_new] = combined_step_project(P_all,F0,V_coarse,F_coarse,delta_t,varargin)
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

  function G = displacement_gradient(V,V0)
    % Gradient of total displacement from (V0,F) to (V,F)
    %
    % Inputs:
    %   V  #V by 3 current positions
    %   V0  #V by 3 positions to take displacement with respect to (e.g. initial
    %     positions)
    G = V0-V;
  end
  function E = displacement_energy(V,V0)
    G = displacement_gradient(V,V0);
    E = trace(G'*G);
  end

  function G = volume_gradient(V,F)
    % Gradient of volume of mesh (V,F)
    G = area_weighted_normal(V,F);
  end
  function E = volume_energy(V,F)
    [~,E] = centroid(V,F);
  end

  function G = symmetry_x_gradient(V,sym_pairs)
    % Gradient of energy compairing differences of pairs of vertices (after
    % reflection due to symmetry over yz-plane)
    %
    % Inputs:
    %   V  #V by 3 current positions
    %   sym_pairs  #sym_pairs by 2 list of indices into V of symmetric pairs
    G = zeros(size(V));
    G(sym_pairs(:,1),1) = 2*(V(sym_pairs(:,1),1)+V(sym_pairs(:,2),1));
    G(sym_pairs(:,1),2) = 2*(V(sym_pairs(:,1),2)-V(sym_pairs(:,2),2));
    G(sym_pairs(:,1),3) = 2*(V(sym_pairs(:,1),3)-V(sym_pairs(:,2),3));
    G(sym_pairs(:,2),1) = grad_sym(sym_pairs(:,2),1) + 2*(V(sym_pairs(:,1),1)+V(sym_pairs(:,2),1));
    G(sym_pairs(:,2),2) = grad_sym(sym_pairs(:,2),2) - 2*(V(sym_pairs(:,1),2)-V(sym_pairs(:,2),2));
    G(sym_pairs(:,2),3) = grad_sym(sym_pairs(:,2),3) - 2*(V(sym_pairs(:,1),3)-V(sym_pairs(:,2),3));
  end
  function E = symmetry_x_energy(V,sym_pairs)
    G = symmetry_x_gradient(V,sym_pairs)
    E = trace(G'*G);
  end

  function G = volumetric_arap_gradient(TV0,TT,TV,delta_t)
    % Gradient of volumetric arap
    % 
    % Inputs:
    %   TV0  #TV>#V by 3 list of tet mesh vertices in initial positions: surface
    %     vertices come first TV0 = [V0;steiner]
    %   TT  #TT by 4 list of tet mesh indices into TV
    %   TV  #TV by 3  new tet mesh vertices (with surface coming first
    %   delta_t  time step used for "conservation of momentum" term

    % Compute one small step of arap flow
    % TV_flow = arap(...)
    G = (TV-TV_flow)/delta_t;
    % Double check sign and magnitude 
  end
  function E = volumetric_arap_energy(TV0,TT,TV)
    % compute current energy of (TV,TT) w.r.t. (TV0,TT)
  end


  disp('started combined_step_project');
  etienne_called = 0;

  simulation_steps = 1;
  energy = 'displacement_initial_and_volume';
  min_progress = 5e-4;
  method = 'shrink_fine_and_expand';
  % define target cage (many times initial mesh)
  CV_0 = V_coarse;
  eps_distance = 1e-4;
  beta_init = 10;
  % parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'simulation_steps'
              ii = ii+1;
              assert(ii<=numel(varargin));
              simulation_steps = varargin{ii};
          case 'energy'
              ii = ii+1;
              assert(ii<=numel(varargin));
              energy = varargin{ii};
          case 'min_progress'
              ii = ii+1;
              assert(ii<=numel(varargin));
              min_progress = varargin{ii};
          case 'method'
              ii = ii+1;
              assert(ii<=numel(varargin));
              method = varargin{ii};
          case 'CV_target'
              ii = ii+1;
              assert(ii<=numel(varargin));
              CV_0 = varargin{ii};
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
  % start with no faces to refine (only needed for 
  % method='shrink_coarse_refining')
  F_to_refine = [];

  % activates the camera toolbar in the figure and set principal rotation
  % axis to none
  cameratoolbar;
  cameratoolbar('SetCoordSys','none');

  axis equal;
  pc = [];
  pv = [];
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])

  % input('Ready. Set. ');
  % fprintf('Go!\n');

  cla;

  % number of steps
  k = size(P_all,3);
  % faces for all vertices of both meshes (for collision detection)
  F_all = [F;CF+size(V,1)];

  % Volumetric ARAP not working, re-visit if needed
  % if strcmp(energy,'volumetric_arap')
  %     % if (volumetric_arap), tetrahedralize initial mesh (rest pose)
  %     [SV_0,SF_0] = tetgen(CV_0,CF);
  % end
  % minimization parameters
  beta = beta_init;
  tol = min_progress;

  % Eltopo parameters
  eps_proximity = eps_distance;
  tol_dt = 1e-1;
  max_min_attempt = 3;
  V_eltopo = [V;CV];
  V_all_prev = V_eltopo;

  if strcmp(energy,'symmetry_x')
      % brute-force way of detecting pairs of points symmetric w.r.t.
      % to x=0 plane
      sym_pairs = [];
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
  CV_filtered = CV;
  CV_orig = CV;
  % starting to time expansion reversing flow: we'll updated CV_filtered after
  % each reverse flow time step.
  tic;
  for t = k-1:-1:1
    % Previous positions of fine mesh
    V_prev = P_all(:,:,t+1);
    % Desired positions of fine mesh after step
    V = P_all(:,:,t);

      switch energy
      case {'displacement_step','displacement_step_and_volume'}
        energy_gradient = @(CV_prev) displacement_gradient(CV_prev,CV_filtered);
        energy_value    = @(CV_prev)   displacement_energy(CV_prev,CV_filtered);
      case {'displacement_initial','displacement_initial_and_volume'}
        energy_gradient = @(CV_prev) displacement_gradient(CV_prev,CV_orig);
        energy_value    = @(CV_prev)   displacement_energy(CV_prev,CV_orig);
      case {'volume'}
        energy_gradient = @(CV_prev) volume_gradient(CV_prev,CF);
        energy_value    = @(CV_prev)   volume_energy(CV_prev,CF);
      case {'symmetry_x'}
        energy_gradient = @(CV_prev) symmetry_x_gradient(CV_prev,sym_pairs);
        energy_value    = @(CV_prev)   symmetry_x_energy(CV_prev,sym_pairs);
      case {'volumetric_arap'}
        % Optimize for internal steiner points according to volumetric arap,
        % fixing surface to CV
        %
        % Note: this is done only once at beginning of this flow step. It is
        % **not** done for each call to `energy_gradient`
        %
        TV = arap(TV0,TT,1:size(CV,1),CV_filtered);
        TV_steiner = TV(size(CV,1)+1:end,:);
        energy_gradient = @(CV_prev) volumetric_arap_gradient(TV0,TT,[CV_prev;TV_steiner],delta_t);
        energy_value  = @(CV_prev)     volumetric_arap_energy(TV0,TT,[CV_prev;TV_steiner]);
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


      %% TODO: THIS SHOULD BE A SEPARATE FUNCTION CALL SO THAT FINAL OPTIMIZATION
      %CAN ALSO USE IT

      E_val = inf;
      % Stepping in energy gradient direction until converged
      bb_iter = 1;
      beta_orig = beta;
      BETA_MIN = 1e-1;
      D_CV_MIN = 1e-4;
      CV_prev = CV_filtered;
      while true
        % Update gradient on coarse mesh
        [CV_grad] = energy_gradient(CV_prev);

        % Try to take one step using el topo
        assert(isempty(intersect_other(V_prev,F,CV_prev,CF,'FirstOnly',true)));
        [V_all_bb,rest_dt] = collide_eltopo_mex( ...
          [V_prev;CV_prev             ],F_all, ...
          [     V;CV_prev-beta*CV_grad], ...
          size(V,1),10*eps_proximity,tol_dt);
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

        % Stop if the change in positions is tiny
        d_CV = max(normrow(CV_filtered - CV_prev));
        if d_CV < D_CV_MIN && bb_iter > 1
          fprintf('Max change in CV (%g) less than D_CV_MIN (%g)\n', ...
            d_CV,D_CV_MIN);
          break;
        end 

        % Compute energy at filtered positions
        E_val_prev = E_val;
        E_val = energy_value(CV_filtered);
        % Is energy decreasing (and not first run)
        if E_val < E_val_prev
          if bb_iter > 1
            % try to increase beta
            beta = min(1.1*beta,beta_orig);
          end
          fprintf('Progress!\n');
        else
          assert(bb_iter > 1);
          % otherwise decrease beta
          beta = beta*0.5;
          % and roll back 
          CV_filtered = CV_prev;
          fprintf('No progress...\n');
        end
        assert(isempty(intersect_other(V,F,CV_filtered,CF,'FirstOnly',true)));

        % Stop if beta is now too small
        if beta < BETA_MIN
          fprintf('Beta (%g) less than BETA_MIN (%g)\n',beta,BETA_MIN);
          break;
        end

        if debug
          tsurf(F,V_prev,'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','b');
          hold on;
          tsurf(CF,CV_filtered,'FaceAlpha',0.5,'EdgeAlpha',0.5,'CData', ...
            sum(abs(CV_prev-CV_filtered).^2,2));
          hold off;
          colorbar
          axis equal;
          title(sprintf('t: %d',t),'FontSize',30);
          drawnow;
        end

        % Updated previous positions for next iteration of loop
        CV_prev = CV_filtered;

        bb_iter = bb_iter+1;
      end
      assert(bb_iter >= 1,'must take at least one step');
      beta = beta_orig;
  end
  time_expansion = toc;

  % TODO: ADD BACK FINAL OPTIMIZATION PHASE?
  time_final_energy = 0;

  V_new = [V;CV_filtered];
  V_coarse_final = CV_filtered;
end 
