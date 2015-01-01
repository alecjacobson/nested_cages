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
  
  quadrature_order = 2;
  V_coarse = [];
  F_coarse = [];
  Pall = [];
  % below parameter only used for ElTopo
  beta_init = 1e-1;
  
  % save timings
  timing.decimation = 0.0;
  timing.flow = 0.0;
  timing.simulation = 0.0;
  
    % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'quadrature_order'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              quadrature_order = varargin{ii};
          case 'beta_init'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              beta_init = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % number of levels
  num_levels = size(levels,2);
  
  % Define last level as input mesh
  V_coarse{2*num_levels+1} = V0;
  F_coarse{2*num_levels+1} = F0;
  
  cages_V{2*num_levels+1} = V0;
  cages_F{2*num_levels+1} = F0;
  
  % loop over different levels
  for k=num_levels:-1:1
      
      tic
      [V_coarse{2*k},F_coarse{2*k}] = cgal_simplification(cages_V{2*k+1},cages_F{2*k+1},levels(k));
      [V_coarse{2*k},F_coarse{2*k}] = meshfix(V_coarse{2*k},F_coarse{2*k});
      timing.decimation = timing.decimation + toc;
      
      
      cla;
      % shirnk fine mesh
      tic
      [Pall,~,F_exp,~] = shrink_fine_expand_coarse_3D(cages_V{2*k+1},cages_F{2*k+1},...
          V_coarse{2*k},F_coarse{2*k},'quadrature_order',quadrature_order,'step_size',5e-1);
      timing.flow = timing.flow + toc;
      Pall_all_times{2*k} = Pall;
      
      % push coarse mesh with physical simulation to obtain the cages
      tic
      % first expand
      pc = [];
      pv = [];
      V = Pall(:,:,end);
      F = cages_F{2*k+1};
      CV = V_coarse{2*k};
      CF = F_coarse{2*k};
      F_all = [F;size(V,1)+CF];
      V_eltopo = [V;CV];
      for j=1:size(Pall,3)-1
          V_all_prev = V_eltopo;
          V = Pall(:,:,end-j);
          fprintf('number of known vertices: %d\n',size(V,1));
          disp('Eltopo expanding with 2*eps')
          [V_eltopo,rest_dt] = collide_eltopo_mex(V_all_prev,F_all,[V;CV],size(V,1),2*sep_thick,1e-1);
          if (rest_dt>0.0)
              disp('ElTopo could not handle it, switching to velocityfilter')
              V_all_prev = inflate_mex(V_all_prev,F_all,size(V,1),eps_proximity);
              V_eltopo = velocity_filter_mex(V_all_prev,[V;CV],F_all,size(V,1),sep_thick,0.01*sep_thick);
          end
          
          CV = V_eltopo(size(V,1)+1:end,:);
          % plot partial results
          cla;
          hold on;
          pc = trisurf(CF,CV(1:end,1),CV(1:end,2),CV(1:end,3),'FaceColor',[0.0 0.0 0.0],'FaceAlpha',0.1);
          pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
          title(sprintf(' flow step: %d/%d', ...
              j, size(Pall,3)-1));
          drawnow;
      end
      fprintf('max difference between simulation and constrained mesh: %g %g %g\n',max(abs(V_eltopo(1:size(V,1),:)-V)));
      disp('done expanding')
      
      % if separation was not achieved, inflate mesh
      V_eltopo = inflate_mex(V_eltopo,F_all,size(V,1),sep_thick);
%       fprintf('distance between meshes before minimization = %g\n', self_distance_mex(V_eltopo,F_all,size(V,1)));
      disp('done inflating')
      
      % now volume minimization
      beta = beta_init;
      
      % configuration is already feasible
      CV = V_eltopo(size(V,1)+1:end,:);
      [~,min_energy] = centroid(CV,CF);
      CV_opt = CV;
      fprintf('initial energy = %g\n', min_energy);
      
      while(beta>0.01*beta_init)
%       while(false)
          
          V_all_prev = V_eltopo;
          
%           areas = doublearea(CV,CF)/2;
%           N = normalizerow(normals(CV,CF));
%           grad_vol = zeros(size(CV));
%           for i=1:size(CV,1)
%               face_idx = mod(find(CF==i)-1,size(CF,1))+1;
%               grad_vol(i,:) = sum([areas(face_idx) areas(face_idx) areas(face_idx)].*N(face_idx,:));
%           end
          grad_vol = area_weighted_normal(CV,CF);
          
          fprintf('number of known vertices: %d\n',size(V,1));
          %                 fprintf('distance between meshes before simulation = %g\n', self_distance_mex(V_all_prev,F_all,size(V,1)));
          %                 input('')
          % step+project
          [Z,rest_dt] = collide_eltopo_mex(V_all_prev,F_all,[V;CV-beta*grad_vol],size(V,1),sep_thick,1e-1);
          if (rest_dt>0.0)
              disp('ElTopo could not handle it, switching to velocityfilter')
              if (~inflated)
                  V_all_prev = inflate_mex(V_all_prev,F_all,size(V,1),sep_thick);
                  % also has to re-calculate the step towards
                  % gradient direction
                  CV_inf = V_all_prev(size(V,1)+1:end,:);
                  areas = doublearea(CV_inf,CF)/2;
                  N = normalizerow(normals(CV_inf,CF));
                  grad_vol = zeros(size(CV_inf));
                  for i=1:size(CV_inf,1)
                      face_idx = mod(find(CF==i)-1,size(CF,1))+1;
                      grad_vol(i,:) = sum([areas(face_idx) areas(face_idx) areas(face_idx)].*N(face_idx,:));
                  end
                  inflated = 1;
              end
              disp('Running velocity_filter_mex');
              Z = velocity_filter_mex(V_all_prev,[V;CV_inf-beta*grad_vol],F_all,size(V,1),sep_thick,0.01*sep_thick);
          end
          
          % if separation was not achieved, inflate mesh
          Z = inflate_mex(Z,F_all,size(V,1),sep_thick);
%           fprintf('distance between meshes before minimization = %g\n', self_distance_mex(Z,F_all,size(V,1)));
          disp('done inflating')
          
          gamma = 1;
          CV_attempt = CV+gamma*(Z(size(V,1)+1:end,:)-CV);
          
          [~,cur_energy] = centroid(CV_attempt,CF);
          
          % if energy decreased, 'save' this state
          if (cur_energy<min_energy)
              % update meshes
              V_eltopo = Z;
              CV = CV_attempt;
              % save optimal result to output
              CV_opt = CV;
              min_energy = cur_energy;
              beta = beta_init;
              inflated = 0;
              fprintf('energy decreased: energy = %g\n', cur_energy);
          else
              beta = 0.5*beta;
              fprintf('energy increased to energy = %g. Changing beta = %g\n', cur_energy, beta);
          end
          % plot partial results
          cla;
          hold on;
          pc = trisurf(CF,CV(1:end,1),CV(1:end,2),CV(1:end,3),'FaceColor',[0.0 0.0 0.0],'FaceAlpha',0.1);
          pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
          title(sprintf('Volume minimization'));
          
          drawnow;
                    
      end
      timing.simulation = timing.simulation + toc;
      
      % output level
      cages_F{2*k} = F_exp;
      cages_V{2*k} = CV_opt;
      V_coarse{2*k} = cages_V{2*k};
      F_coarse{2*k} = cages_F{2*k};
      
      % save partial
      save('matrioshka_partial.mat','Pall','cages_V','cages_F');
      
      fprintf('done generated layer. Now generation of the wall');
      % first copy layer
      V_coarse{2*k-1} = V_coarse{2*k};
      F_coarse{2*k-1} = F_coarse{2*k};
      % shirnk fine mesh
      cla;
      tic
      [Pall,~,F_exp,~] = shrink_fine_expand_coarse_3D(cages_V{2*k},cages_F{2*k},...
          V_coarse{2*k-1}+1e-6*randn(size(V_coarse{2*k-1})),F_coarse{2*k-1},'quadrature_order',quadrature_order,'step_size',5e-1);
      timing.flow = timing.flow + toc;
      Pall_all_times{2*k-1} = Pall;
      timing.flow = timing.flow + toc;
      
      % push coarse mesh with physical simulation to obtain the cages
      tic
      % first expand
      pc = [];
      pv = [];
      V = Pall(:,:,end);
      F = cages_F{2*k};
      CV = V_coarse{2*k-1};
      CF = F_coarse{2*k-1};
      F_all = [F;size(V,1)+CF];
      V_eltopo = [V;CV];
      for j=1:size(Pall,3)-1
          V_all_prev = V_eltopo;
          V = Pall(:,:,end-j);
          fprintf('number of known vertices: %d\n',size(V,1));
          disp('Eltopo expanding with 2*eps')
          [V_eltopo,rest_dt] = collide_eltopo_mex(V_all_prev,F_all,[V;CV],size(V,1),2*wall_thick,1e-1);
          if (rest_dt>0.0)
              disp('ElTopo could not handle it, switching to velocityfilter')
              V_all_prev = inflate_mex(V_all_prev,F_all,size(V,1),wall_thick);
              V_eltopo = velocity_filter_mex(V_all_prev,[V;CV],F_all,size(V,1),wall_thick,0.01*wall_thick);
          end
          
          CV = V_eltopo(size(V,1)+1:end,:);
          % plot partial results
          cla;
          hold on;
          pc = trisurf(CF,CV(1:end,1),CV(1:end,2),CV(1:end,3),'FaceColor',[0.0 0.0 0.0],'FaceAlpha',0.1);
          pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
          title(sprintf(' flow step: %d/%d', ...
              j, size(Pall,3)-1));
          drawnow;
      end
      fprintf('max difference between simulation and constrained mesh: %g %g %g\n',max(abs(V_eltopo(1:size(V,1),:)-V)));
      disp('done expanding')
      
      % if separation was not achieved, inflate mesh
      V_eltopo = inflate_mex(V_eltopo,F_all,size(V,1),wall_thick);
%       fprintf('distance between meshes before minimization = %g\n', self_distance_mex(V_eltopo,F_all,size(V,1)));
      disp('done inflating')
      
      % now volume minimization
      beta = beta_init;
      
      % configuration is already feasible
      CV = V_eltopo(size(V,1)+1:end,:);
      [~,min_energy] = centroid(CV,CF);
      CV_opt = CV;
      fprintf('initial energy = %g\n', min_energy);
      
      while(beta>0.01*beta_init)
%       while(false)
          
          V_all_prev = V_eltopo;
          
%           areas = doublearea(CV,CF)/2;
%           N = normalizerow(normals(CV,CF));
%           grad_vol = zeros(size(CV));
%           for i=1:size(CV,1)
%               face_idx = mod(find(CF==i)-1,size(CF,1))+1;
%               grad_vol(i,:) = sum([areas(face_idx) areas(face_idx) areas(face_idx)].*N(face_idx,:));
%           end
          grad_vol = area_weighted_normal(CV,CF);
          
          fprintf('number of known vertices: %d\n',size(V,1));
          %                 fprintf('distance between meshes before simulation = %g\n', self_distance_mex(V_all_prev,F_all,size(V,1)));
          %                 input('')
          % step+project
          [Z,rest_dt] = collide_eltopo_mex(V_all_prev,F_all,[V;CV-beta*grad_vol],size(V,1),wall_thick,1e-1);
          if (rest_dt>0.0)
              disp('ElTopo could not handle it, switching to velocityfilter')
              if (~inflated)
                  V_all_prev = inflate_mex(V_all_prev,F_all,size(V,1),wall_thick);
                  % also has to re-calculate the step towards
                  % gradient direction
                  CV_inf = V_all_prev(size(V,1)+1:end,:);
                  areas = doublearea(CV_inf,CF)/2;
                  N = normalizerow(normals(CV_inf,CF));
                  grad_vol = zeros(size(CV_inf));
                  for i=1:size(CV_inf,1)
                      face_idx = mod(find(CF==i)-1,size(CF,1))+1;
                      grad_vol(i,:) = sum([areas(face_idx) areas(face_idx) areas(face_idx)].*N(face_idx,:));
                  end
                  inflated = 1;
              end
              disp('Running velocity_filter_mex');
              Z = velocity_filter_mex(V_all_prev,[V;CV_inf-beta*grad_vol],F_all,size(V,1),wall_thick,0.01*wall_thick);
          end
          
          % if separation was not achieved, inflate mesh
          Z = inflate_mex(Z,F_all,size(V,1),wall_thick);
%           fprintf('distance between meshes before minimization = %g\n', self_distance_mex(Z,F_all,size(V,1)));
          disp('done inflating')
          
          gamma = 1;
          CV_attempt = CV+gamma*(Z(size(V,1)+1:end,:)-CV);
          
          [~,cur_energy] = centroid(CV_attempt,CF);
          
          % if energy decreased, 'save' this state
          if (cur_energy<min_energy)
              % update meshes
              V_eltopo = Z;
              CV = CV_attempt;
              % save optimal result to output
              CV_opt = CV;
              min_energy = cur_energy;
              beta = beta_init;
              inflated = 0;
              fprintf('energy decreased: energy = %g\n', cur_energy);
          else
              beta = 0.5*beta;
              fprintf('energy increased to energy = %g. Changing beta = %g\n', cur_energy, beta);
          end
          
          % plot partial results
          hold on;
          cla;
          pc = trisurf(CF,CV(1:end,1),CV(1:end,2),CV(1:end,3),'FaceColor',[0.0 0.0 0.0],'FaceAlpha',0.1);
          pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
          title(sprintf('Volume minimization'));
          
          drawnow;
                    
      end
      timing.simulation = timing.simulation + toc;
      
      % output level
      cages_F{2*k-1} = F_exp;
      cages_V{2*k-1} = CV_opt;
      V_coarse{2*k-1} = cages_V{2*k-1};
      F_coarse{2*k-1} = cages_F{2*k-1};
      
      % save partial
      save('matrioshka_partial.mat','Pall','cages_V','cages_F');
            
      
  end
  
end