function [V_disp,dt,violations,IF] = signed_distance_normal_displacement(V0,F0,varargin)
%
% Given an initial mesh (V0,F0), computes normal displacements
% V0 = V0+dt*(-N) until it finds a normal displacement that does
% not intersect the original. If it doesn't, output dt = +inf.
%
% [V_disp,dt,violations] = signed_distance_normal_displacement(V0,F0)
%
% Input:
% (V0,F0): Input mesh
% Output
% V_disp:    = V0+dt*(-N) such that (V_disp,F0) doesn't intersect (V0,F0)
% dt:     dt above
% violations:   number of times that a normal at a vertex has dot
% product with some face around it smaller than zero

ii = 1;
vis = false;
while ii < numel(varargin)
    switch varargin{ii}
        case 'vis'
            ii = ii+1;
            assert(ii<=numel(varargin));
            vis = varargin{ii};
        otherwise
            error('Unsupported parameter: %s',varargin{ii});
    end
    ii = ii+1;
end

% calculate normals according to signed_distance_direction function
[N_dist,~,~] = signed_distance_direction(V0,V0,F0);

% Triangle normals
N_faces = normals(V0,F0);

if (vis)
    cameratoolbar;
    
    % plot general options
    axis equal;
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    % initialize plot handles
    % clear title
    title('shrinking','FontSize',30);
    hold on;
    
    pc = trisurf(F0,V0(:,1),V0(:,2),V0(:,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.05,'EdgeAlpha',0.2);
    pv = trisurf(F0,V0(:,1),V0(:,2),V0(:,3),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.4,'EdgeAlpha',0.2);
    % draw everything
    plot_struct = struct('V_shrink',V0,'F_shrink',F0,'V_to_intersect',V0,...
      'F_to_intersect',F0,'V_coarse',V0,'F_coarse',F0, 'quad_points',[],...
      'grad_vertices',[],'grad_quad',[],'quad_closest',[],'V_shrink_prev',V0,'IF',[],'pc',pc,'pc1',[],'pc_alpha',0.05,'pv',pv,'pv_alpha',0.4,...
      'p_quiver',[],'p_quiver_quad',[],'p_close',[],'p_quad',[],'show_weights',0);
    
%     plot_strict.p_quiver = quiver3(V0(:,1),V0(:,2),V0(:,3),...
%           0.1*N_dist(:,1),0.1*N_dist(:,2),0.1*N_dist(:,3),0);
%       set(plot_struct.p_quiver,'Color',[0 0.75 0.0],'Linewidth',1.5);
      
    drawnow;
  
end

% compute number of violations
violations = 0;
for i=1:size(V0,1)
    faces = mod(find(F0 == i)-1,size(F0,1))+1;
    for k=1:size(faces,1)
        if (dot(N_faces(faces(k),:),N_dist(i,:))>0)
            violations = violations+1;
        end
    end
end

% Now compute normal displacements
dt = 1e-32;
IF = 1;
while (size(IF,1)>0 && dt<=1e-2)
    
%     V_disp = V0;
    
%     while (any(V_disp(:)==V0(:)))
%         dt = dt*2
%         V_disp = V0 + dt*N_dist;
%     end
    
    V_disp = V0 + dt*N_dist;
    if (~any(V_disp(:)==V0(:)))
        IF = intersect_other(V0,F0,V_disp,F0);
%         input('');
        if (vis)
            plot_struct.V_shrink = V_disp;
            plot_struct.F_shrink = F0;
            plot_struct = flow_plot_control(plot_struct,true);
        end
    end
    size(IF)
    IF
    dt = dt*1.2
%     meshplot(V_disp,F0([361 552],:))
    
end
% If couldn't find dt, output infinity
if (dt>1e-2)
    V_disp = zeros(size(V0));
    dt = inf;
end

