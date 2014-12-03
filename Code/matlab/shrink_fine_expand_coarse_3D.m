function [Pall_fine,Pall_coarse,F_exp,F_shrink] = shrink_fine_expand_coarse_3D(V0,F0,V_coarse,F_coarse,varargin)
  % SHRINK_FINE_EXPAND_COARSE_3D
  % [Pall_fine,Pall_coarse] = shrink_fine_expand_coarse_3D(V0,F0,V_coarse,F_coarse,varargin)
  %
  % Add description
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
  %   coarse mesh
  %   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
  %   of the coarse mesh
  %   Optional:
  %     'quadrature_order': 1, 2 or 3 (default=3)
  % Output:
  %   Pall_fine   (#vertices)x3xsteps lists of mesh vertex positions of 
  %   the fine mesh at each time step
  %   Pall_coarse   (#vertices_cage)x3xsteps lists of mesh vertex positions of 
  %   the fine mesh at each time step
  
%   % refine faces connected to problematic vertices
%   face_idx = union(mod(find(F_coarse==209)-1,size(F_coarse,1))+1,mod(find(F_coarse==197)-1,size(F_coarse,1))+1);
%   [V_coarse,F_coarse] = upsample(V_coarse,F_coarse,'OnlySelected',face_idx);
  

  % plot general options
  cla;
  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])
  hold on
  cameratoolbar;
  cameratoolbar('SetCoordSys','none');

% %   attempt 3: consider gradients coming from quadrature calculations
% 
%     % calculate quadrature-based gradients using during flow
%     % the call below does much more things then needed, do a separate 
%     % simpler function if needed
%     for b=1:1
%         [~,~,grad_coarse,~,~,~,~] = signed_distance_direction_quadrature_matrix(V_coarse,F_coarse,1,3,V_coarse,F_coarse,...
%             1:size(V_coarse,1),'scalar_search',0.5*doublearea(V_coarse,F_coarse),...
%             grad_quadrature_to_vertices(V_coarse,F_coarse,0.5*doublearea(V_coarse,F_coarse),3),massmatrix(V_coarse,F_coarse,'barycentric'),...
%             0.0,cotmatrix(V_coarse,F_coarse),1);
%         
% %         p_coarse = trisurf(F_coarse,V_coarse(:,1),V_coarse(:,2),V_coarse(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.05);
% %         p_quiver_1 = quiver3(V_coarse(:,1),V_coarse(:,2),V_coarse(:,3),100*grad_coarse(:,1),100*grad_coarse(:,2),100*grad_coarse(:,3),0);
% %         set(p_quiver_1,'Color','m','Linewidth',1.5);
% %         input('');
% %         delete(p_quiver_1);
% %         delete(p_coarse);
%         prob_vertices = [];
%         for a=1:size(V_coarse,1)
%             if norm(grad_coarse(a,:))<5e-4
%                 prob_vertices = [prob_vertices a];
% %                 a
% %                 p_quiver_2 = quiver3(V_coarse(a,1),V_coarse(a,2),V_coarse(a,3),100*grad_coarse(a,1),100*grad_coarse(a,2),100*grad_coarse(a,3),0);
% %                 set(p_quiver_2,'Color','m','Linewidth',1.5);
% %                 input('');
% %                 delete(p_quiver_2);
%             end
%         end
%     %     prob_vertices
%     %     input('');
%         % refine faces connected to problematic vertices
%         face_idx = [];
%         for a=1:size(prob_vertices,2)
%             face_idx = union(mod(find(F_coarse==prob_vertices(a))-1,size(F_coarse,1))+1,face_idx);
%         end
%     %     face_idx
%     %     input('');
%         if (~isempty(face_idx))
%             [V_coarse,F_coarse] = upsample(V_coarse,F_coarse,'OnlySelected',face_idx);
%         end
%     end
  

  
  
  quadrature_order = 3;
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'quadrature_order'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              quadrature_order = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % initialize (V_exp,F_exp) as (V_coarse,F_coarse)
  V_exp = V_coarse;
  F_exp = F_coarse;
  
%   tsurf(F_exp,V_exp,'FaceIndices',true,'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.5);
  
  % initialize (V_shrink,F_shrink) as (V0,F0)
  V_shrink = V0;
  F_shrink = F0;
  
%   % upsample twice
%   [V_shrink,F_shrink] = upsample(V_shrink,F_shrink,'OnlySelected',1:size(F_shrink,1));
%   [V_shrink,F_shrink] = upsample(V_shrink,F_shrink,'OnlySelected',1:size(F_shrink,1));
  
  % concatenate meshes and test for intersections
  V_all = [V_exp;V_shrink];
  F_all = [F_exp;F_shrink+size(V_exp,1)];
  [~,~,IF] = selfintersect(V_all,F_all,'DetectOnly',true);
  
  %initialize flow
  ii = 1;
  Pall_coarse(:,:,1) = V_exp;
  Pall_fine(:,:,1) = V_shrink;
  
  
  % flow settings for the coarse mesh
  area_initial_coarse = doublearea(V_coarse,F_coarse)/2;
  A_qv_coarse = grad_quadrature_to_vertices(V_coarse,F_coarse,area_initial_coarse,quadrature_order);
  % needed to run 'signed_distance_direction_quadrature_matrix'
  % remove in the future
  M_coarse = massmatrix(V_coarse,F_coarse,'barycentric');
  w_lap_coarse = 0.0;
  L_coarse = cotmatrix(V_coarse,F_coarse);
  % initialize all vertices as moving vertices
  moving_vertices_coarse = ones(size(V_exp,1),1);
  
  % flow settings for the fine mesh
  area_initial_shrink = doublearea(V_shrink,F_shrink)/2;
  A_qv_shrink = grad_quadrature_to_vertices(V_shrink,F_shrink,area_initial_shrink,quadrature_order);
  % needed to run 'signed_distance_direction_quadrature_matrix'
  % remove in the future
  M_shrink = massmatrix(V_shrink,F_shrink,'barycentric');
  w_lap_shrink = 0.0;
  L_shrink = cotmatrix(V_shrink,F_shrink);
  % initialize all vertices as moving vertices
  moving_vertices_shrink = ones(size(V_shrink,1),1);
  
  % plot general options
  axis equal;
  set(gca,'xtick',[])
  set(gca,'xticklabel',[])
  set(gca,'ytick',[])
  set(gca,'yticklabel',[])
  
  % initialize plot handles
  hold on;
  pc = trisurf(F_exp,V_exp(:,1),V_exp(:,2),V_exp(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.05);
  pv = trisurf(F_shrink,V_shrink(:,1),V_shrink(:,2),V_shrink(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.4);
  % draw everything
  drawnow;
  
  % define plotting struct
  plot_struct = struct('V_shrink',V_shrink,'F_shrink',F_shrink,'V_to_intersect',V_exp,...
      'F_to_intersect',F_exp,'V_coarse',V_exp,'F_coarse',F_exp, 'quad_points',[],...
      'grad_vertices',[],'grad_quad',[],'quad_closest',[],'V_shrink_prev',V0,'IF',[],'pc',pc,'pc1',[],'pc_alpha',0.05,'pv',pv,'pv_alpha',0.4,...
      'p_quiver',[],'p_quiver_quad',[],'p_close',[],'p_quad',[],'show_weights',0);
  
  % this fixes the axis for the plots
  axis manual;
  
  % activates the camera toolbar in the figure and set principal rotation
  % axis to none
  cameratoolbar;
  cameratoolbar('SetCoordSys','none');
  
  %   input('');
  hold off;
  flow_attempt = 1;
  
  while(size(IF,1)>0 || flow_attempt==1)
      
      if (flow_attempt>1)
          F_shrink_prev = F_shrink;
          [V_shrink,F_shrink] = upsample(Pall_fine(:,:,end),F_shrink,'OnlySelected',IF(:,1));
          clear Pall_fine_new;
          for b=1:size(Pall_fine,3)-1
              [Pall_fine_new(:,:,b),~] = upsample(Pall_fine(:,:,b),F_shrink_prev,'OnlySelected',IF(:,1));
          end
          clear Pall_fine;
          Pall_fine = Pall_fine_new;
          Pall_fine(:,:,end+1) = V_shrink;
          % flow settings for the fine mesh
          area_initial_shrink = doublearea(V_shrink,F_shrink)/2;
          A_qv_shrink = grad_quadrature_to_vertices(V_shrink,F_shrink,area_initial_shrink,quadrature_order);
          % needed to run 'signed_distance_direction_quadrature_matrix'
          % remove in the future
          M_shrink = massmatrix(V_shrink,F_shrink,'barycentric');
          w_lap_shrink = 0.0;
          L_shrink = cotmatrix(V_shrink,F_shrink);
          % initialize all vertices as moving vertices
          moving_vertices_shrink = ones(size(V_shrink,1),1);
          ii = ii+1;
%           ii = 1;
%           clear Pall_fine;
%           Pall_fine(:,:,1) = V_shrink;
      end
      flow_attempt = flow_attempt+1;
  
      while (size(IF,1)> 0 && mod(ii,100)~=0)      

          ii = ii + 1

          % first one step to expand coarse mesh
          V_exp_prev = V_exp;

%           if (mod(ii,10)==0)
          if (false)
              [V_exp,moving_vertices_coarse,grad_energy,quad,grad_quadrature,C,~] = signed_distance_direction_quadrature_matrix(V_exp,F_exp,1,quadrature_order,...
                  V_coarse,F_coarse,moving_vertices_coarse,'scalar_search',area_initial_coarse,A_qv_coarse,M_coarse,w_lap_coarse,L_coarse,-1);

              % Eltopo to remove intersections on the coarse mesh
              [V_exp,~] = collide_eltopo_mex(V_exp_prev,F_exp,V_exp,0,1e-4,1e-10);

              Pall_coarse(:,:,ii) = V_exp;

              IF = intersect_other(V_shrink,F_shrink,V_exp,F_exp);
              disp('number of intersections after coarse mesh expansion')
              size(IF)

              plot_struct.quad_points = quad;
              plot_struct.grad_vertices = grad_energy;
              plot_struct.grad_quad = grad_quadrature;
              plot_struct.quad_closest = C;
              plot_struct.V_coarse = V_exp;
              plot_struct.F_coarse = F_exp;
              plot_struct.V_to_intersect = V_exp;
              plot_struct.F_to_intersect = F_exp;
              %       plot_struct.V_shrink_prev = V_exp_prev;
              plot_struct = flow_plot_control(plot_struct,false);


          else

              % shrinking step

              % shrink fine mesh 5 steps
              V_shrink_prev = V_shrink;
              [V_shrink,moving_vertices_shrink,grad_energy,quad,grad_quadrature,C,~] = signed_distance_direction_quadrature_matrix(V_shrink,F_shrink,1,quadrature_order,...
                  V_coarse,F_coarse,moving_vertices_shrink,'scalar_search',area_initial_shrink,A_qv_shrink,M_shrink,w_lap_shrink,L_shrink,1);


              %               % refine faces that are too streched and intersect the coarse mesh
              %               if mod(ii,10)==0
              %                   IF_fine = unique(IF(:,2));
              %                   int_edges = [F_shrink(IF_fine,1) F_shrink(IF_fine,2);...
              %                       F_shrink(IF_fine,2) F_shrink(IF_fine,3);...
              %                       F_shrink(IF_fine,3) F_shrink(IF_fine,1)];
              %                   initial_edge_lengths = normrow(Pall_fine(int_edges(:,1),:,1)-Pall_fine(int_edges(:,2),:,1));
              %                   current_edge_lengths = normrow(Pall_fine(int_edges(:,1),:,end-1)-Pall_fine(int_edges(:,2),:,end-1));
              %                   stretched_edges = find((current_edge_lengths./initial_edge_lengths)>1.0);
              %                   IF_fine = unique(IF_fine(mod(stretched_edges-1,size(IF_fine,1))+1));
              %
              %                   plot_struct.V_shrink = V_shrink;
              %                   plot_struct.F_shrink = F_shrink;
              %                   plot_struct = flow_plot_control(plot_struct,true);
              %
              %                   if ~isempty(IF_fine)
              %                       F_shrink_old = F_shrink;
              %                       [V_shrink,F_shrink] = upsample(V_shrink,F_shrink,'OnlySelected',IF_fine);
              %                       Pall_new = zeros(size(V_shrink,1),3,ii);
              %                       for t=1:size(Pall_fine,3)-2
              %                           [Pall_new(:,:,t),~] = upsample(Pall_fine(:,:,t),F_shrink_old,'OnlySelected',IF_fine);
              %                       end
              %                       Pall_fine = Pall_new;
              %                   end
              %
              %
              %                   disp('new number of faces')
              %                   size(F_shrink,1)
              %                   disp('new number of vertices')
              %                   size(V_shrink,1)
              %                   input('')
              %
              %                   V_shrink_old = V_shrink;
              %     %               [V_shrink] = shrink_vertices_individually(V_shrink,F_shrink,V_coarse,F_coarse,size(V_shrink_prev,1)+1:size(V_shrink,1));
              %
              %                   area_initial_shrink = doublearea(Pall_fine(:,:,1),F_shrink)/2;
              %                   A_qv_shrink = grad_quadrature_to_vertices(V_shrink,F_shrink,area_initial_shrink,quadrature_order);
              %                   M_shrink = massmatrix(Pall_fine(:,:,1),F_shrink,'barycentric');
              %                   L_shrink = cotmatrix(Pall_fine(:,:,1),F_shrink);
              %                   moving_vertices_shrink = ones(size(V_shrink,1),1);
              %               end

              plot_struct.quad_points = quad;
              plot_struct.grad_vertices = grad_energy;
              plot_struct.grad_quad = grad_quadrature;
              plot_struct.quad_closest = C;
              plot_struct.V_shrink = V_shrink;
              plot_struct.V_shrink_prev = V_shrink_prev;
              plot_struct.F_shrink = F_shrink;
              plot_struct = flow_plot_control(plot_struct,false);

              Pall_fine(:,:,end+1) = V_shrink;
              
              V_shrink_prev = V_shrink;
              IF = intersect_other(Pall_fine(:,:,end),F_shrink,Pall_coarse(:,:,end),F_coarse);
              disp('number of intersections after fine mesh shrinking')
              size(IF)

          end

          %       hold off;


      end
  
  end