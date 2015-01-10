function [CV_filtered,timing] = ...
  one_step_project( ...
    V_prev,V,F,CV,CF, ...
    energy_gradient,energy_value,cb_data, ...
    varargin)
  % Given (V_prev,F) previous fine mesh, (V,F) new fine mesh and
  % (CV,CF) current coarse mesh, step from (V_prev,F) to (V,F) and
  % obtian energy minimizing (CV,CF).
  % 
  % Inputs:
  %   (V_prev,F)  previous fine mesh
  %   (V,F) next fine mesh
  %   (CV,CF) coarse mesh in the begining of the time step
  %   energy_gradient  function handle computing energy gradient
  %   energy_value  function handle computing energy value
  %   cb_data  callback data for the energy function handles
  %   Optional:
  %     'Tol': tolerance for Eltopo stop trying to cut the time step
  %     'PlotInfo': info struct for plotting
  %     'BetaInit'  Initial value for beta the step size
  %     'Eps'  epsilon used as separation distance inside el topo
  %     'Debug'  whether to show debug info
  % Output:
  %   CV_filtered: new embedding for the coarse mesh
  %   timing statistics for this step

  function flag = is_converged(CV_filtered,CV_prev,beta,E_val)
    % Stop if the change in positions is tiny
    d_CV = max(normrow(CV_filtered - CV_prev));
    if d_CV < D_CV_MIN && bb_iter > 1
      fprintf('\nMax change in CV (%g) less than D_CV_MIN (%g)\n', ...
        d_CV,D_CV_MIN);
      flag = true;
      return;
    end 
    if E_val <= 0
      fprintf('\nEnergy (%g) is <=0\n',E_val);
      flag = true;
      return;
    end
    % Stop if beta is now too small
    if beta < BETA_MIN
      fprintf('\nBeta (%g) less than BETA_MIN (%g)\n',beta,BETA_MIN);
      flag = true;
      return;
    end
    flag = false;
  end

  % define target cage (many times initial mesh)
  eps_distance = 1e-4;
  beta_init = 1e-2;
  tol_dt = 1e-1;
  plot_info = [];
  plot_info.t = 0;
  plot_info.energy = 'unknown';
  debug = true;
  skip_el_topo = false;
  D_CV_MIN = 1e-5;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {  'SkipElTopo','Eps','BetaInit','Tol','PlotInfo','Debug','D_CV_MIN'}, ...
    {'skip_el_topo','eps_distance','beta_init','tol_dt','plot_info','debug','D_CV_MIN'});
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

  % initialize CV_filtered
  CV_filtered = CV;
  timing = [];
  timing.etienne_called = 0;
  pc = [];
  pv = [];

  % V_all_prev  where the meshes were
  % [V;CV]   where the meshes want to go
  % V_eltopo  where black box moved the meshes 
  % Then update CV

  % #Parameters
  %   beta  Magnitude of coarse mesh energy gradient we're attempting
  %   eps_distance  desired separation distance (10*eps for el topo, eps
  %      for etienne)
  %   out_proximity  separation distance determined by etienne's inflation
  %


  % faces for all vertices of both meshes (for collision detection)
  F_all = [F;CF+size(V,1)];

  beta = beta_init;
  % initialize energy as inf
  E_val = inf;
  E_opt = inf;
  % Stepping in energy gradient direction until converged
  bb_iter = 1;
  BETA_MIN = 1e-3;
  %BETA_MIN = 1e-6;
  %D_CV_MIN = 1e-6;
  CV_prev = CV_filtered;
  while true
    % Update gradient on coarse mesh
    [CV_grad,cb_data] = energy_gradient(CV_prev,cb_data);
    %[dbE] = energy_value(CV_prev-beta*CV_grad,cb_data);
    %[cb_data.E dbE]
    %error

    assert(isempty(intersect_other(V_prev,F,CV_prev,CF,'FirstOnly',true)));
    [~,~,siIF] = selfintersect(CV_prev,CF,'DetectOnly',true,'FirstOnly',true);
    assert(isempty(siIF));

    while true


      if max(abs(CV_grad(:)))>1
        cla;
        hold on;
        tsurf(CF,CV_prev, ...
            'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.2,'EdgeAlpha',0.2);
        tsurf(CF,CV_prev-CV_grad, ...
            'FaceColor',[0.0 0.5 0.0],'FaceAlpha',0.2,'EdgeAlpha',0.2);
        hold off;
      end
      assert(max(abs(CV_grad(:)))<1,'Gradient is too big. Aborting...');

      % Try to take one step using el topo
      %assert(isempty(intersect_other(V_prev,F,CV_prev,CF,'FirstOnly',true)));
      if ~skip_el_topo
        [V_all_bb,rest_dt] = collide_eltopo_mex( ...
          [V_prev;CV_prev             ],F_all, ...
          [     V;CV_prev-beta*CV_grad], ...
          size(V,1),eps_distance,tol_dt);
          % Did el topo fail?
          if (rest_dt>0.0)
            disp('ElTopo could not handle it, switching to velocityfilter')
            % Try to inflate current situation to accomodate eps_distance
            [V_all_inf,out_proximity] = inflate_mex( ...
              [V_prev;CV_prev],F_all, ...
              size(V,1),eps_distance);
            V_all_bb = velocity_filter_mex( ...
              V_all_inf, ... % current position (after some possible inflation)
              [V;CV_prev - beta*CV_grad], ... % desired position (after beta time step)
              F_all,size(V,1),out_proximity,0.01*out_proximity);
            timing.etienne_called = timing.etienne_called + 1; 
          end
      else
          V_all_bb = velocity_filter_mex( ...
              [V_prev;CV_prev], ...
              [V;CV_prev-beta*CV_grad], ...
              F_all,size(V,1),eps_distance,0.01*eps_distance);
      end
      % At this point, V_all_bb contains vertex positions for both meshes
      % after the "Black Box" velocity filter (el topo + etienne)

      % Extract new positions for coarse mesh
      CV_filtered = V_all_bb(size(V,1)+1:end,:);
      % (CV_filtered,CF) should not intersect (V,F)
      assert(isempty(intersect_other(V,F,CV_filtered,CF,'FirstOnly',true)));
      [~,~,siIF] = selfintersect(CV_filtered,CF,'DetectOnly',true,'FirstOnly',true);
      assert(isempty(siIF));

      [E_val,cb_data] = energy_value(CV_filtered,cb_data);
      % Is energy decreasing (and not first run)
      if E_val < E_opt
        E_opt = E_val;
        if bb_iter > 1
          % try to increase beta
          beta = min(1.1*beta,beta_init);
        end
        fprintf('-');
        break;
      end

      assert(bb_iter > 1);
      % otherwise decrease beta and continue
      beta = beta*0.5;
      fprintf('+');

      if is_converged(CV_filtered,CV_prev,beta,E_val)
        % Use last state since it had less energy
        CV_prev = CV_filtered;
        break;
      end
    end


    %assert(isempty(intersect_other(V,F,CV_filtered,CF,'FirstOnly',true)));

    if debug
      hold on;
        axis equal;
        % delete previous plot
%         delete(pc);
%         delete(pv);
        cla;
        % trisurf maintains previous axes, while tsuyrf doesn't
        pv = trisurf(F,V(:,1),V(:,2),V(:,3),...
            'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2,'EdgeAlpha',0.2);
        pc = trisurf(CF,CV_filtered(:,1),CV_filtered(:,2),CV_filtered(:,3),...
            'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.1,'EdgeAlpha',0.2);
        title(sprintf('energy: %s, t: %d',plot_info.energy,plot_info.t),'FontSize',20,'Interpreter','none');
        drawnow;
      hold off;
    end

    % Stop if the change in positions is tiny
    if is_converged(CV_filtered,CV_prev,inf,E_val)
      break;
    end

    % Updated previous positions for next iteration of loop
    CV_prev = CV_filtered;
    % Assuming we have succeeded in moving the fine mesh, then V_prev
    % should stay put
    V_prev = V;

    bb_iter = bb_iter+1;
  end
  assert(bb_iter >= 1,'must take at least one step');
end
