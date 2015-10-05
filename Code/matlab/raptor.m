[V0,F] = load_mesh('../../Meshes/Results/raptor_claw/claw_input.obj');
[V1,F] = load_mesh('../../Meshes/Results/raptor_claw/claw_flow_standard.obj');
[Q] = load_mesh('../../Meshes/Results/raptor_claw/claw_flow_standard_quad_points.obj');
[CV,CF] = load_mesh('../../Meshes/Results/raptor_claw/claw_coarse_initial.obj');
R = axisangle2matrix([0 1 0],-pi*0.6);
V0 = V0*R;
V1 = V1*R;
Q = Q*R;
CV = CV*R;

BB = [ min([CV;V0]);max([CV;V0])];
data = render_in_cage(V0,F,CV,CF, ...
  'ColorIntersections',false, ...
  'View',[-54,10]);

T = get(gca,'tightinset');set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

drawnow;
imwrite(myaa('raw'),'raptor-claw-00.png');

data = render_in_cage(V1,F,CV,CF,'Data',data,'View',[-54,10]);
data.t.FaceAlpha = 1;

imwrite(myaa('raw'),'raptor-claw-01.png');

%delete(t);
%delete(s);
%
%hold on;
%s = scatter3(Q(:,1),Q(:,2),Q(:,3),'.','SizeData',100);
%hold off;
%ss = copyobj(s,s.Parent);
%light_pos = [l.Position strcmp(l.Style,'local')];
%ground = [0 0 -1 min(NV(:,3))];
%d = ground * light_pos';
%shadow_mat = d*eye(4) - light_pos'*ground;
%U = [Q ones(size(Q,1),1)]*shadow_mat';
%U = bsxfun(@rdivide,U(:,1:3),U(:,4));
%set(ss,'XData',U(:,1),'YData',U(:,2),'ZData',U(:,3));
%ss.MarkerEdgeColor = [0.8 0.8 0.8];
