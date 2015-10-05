addpath ../../Meshes/Results/point_clouds/giraffe-2k-volume/
[V,F] = load_mesh('giraffe_0.obj');
[CV,CF] = load_mesh('giraffe_1.obj');
V = V*axisangle2matrix([1 0 0],-pi/2);
CV = CV*axisangle2matrix([1 0 0],-pi/2);
s = scatter3(V(:,1),V(:,2),V(:,3),'.');
hold on;
tc = tsurf(CF,CV,'FaceAlpha',0.3,'EdgeColor','k','FaceColor',[0.9 0.9 0.9]); view(2); axis equal;
tc.LineWidth = 4;
hold off;

l = light('Position',[-0.1 -0.9 0.6],'Style','infinite');
sc = add_shadow(tc,l,'Ground',[0 0 -1 min(CV(:,3))-1e-4]);
sc.FaceColor = [0.95 0.95 0.95];

ss = copyobj(s,s.Parent);
light_pos = [l.Position strcmp(l.Style,'local')];
ground = [0 0 -1 min(CV(:,3))];
d = ground * light_pos';
shadow_mat = d*eye(4) - light_pos'*ground;
U = [V ones(size(V,1),1)]*shadow_mat';
U = bsxfun(@rdivide,U(:,1:3),U(:,4));
set(ss,'XData',U(:,1),'YData',U(:,2),'ZData',U(:,3));
ss.MarkerEdgeColor = [0.8 0.8 0.8];

camproj('persp');
view(-47,14);
axis equal;
axis manual;
T = get(gca,'tightinset');
set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);

pause
set(gcf,'Color','w');
set(gca,'Visible','off');
imwrite(myaa('raw'),'giraffe-points-with-cage.png');

delete(tc);
delete(sc);

imwrite(myaa('raw'),'giraffe-points.png');
