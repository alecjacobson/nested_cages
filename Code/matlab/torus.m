[V,F] = load_mesh('squiggly-torus.ply');
t = tsurf(F,V*axisangle2matrix([1 0 0],-pi/2), ...
  'FaceAlpha',0.8,'EdgeColor','none', ...
  'FaceVertexCData',repmat([0.0329 0.3724 0.8765],size(V,1),1),fphong);


%gray = [0.840 0.82 0.8];
%[CV,CF] = load_mesh('Csaszar.obj');
%CV(7,:) = [0 0 7.525];
%CV(:,3) = CV(:,3)*0.8;
%t = tsurf(CF, ...
%  CV*axisangle2matrix([1 0 0],-pi/2)*axisangle2matrix([0 0 1],pi*0.48)*axisangle2matrix([1 0 0],-0.3) ...
%  ,'FaceAlpha',0.8,'LineWidth',4,'EdgeColor','k','FaceColor',gray);


l = light('Position',[ 0.8 -0.9 0.8],'Style','infinite');
%s = add_shadow(t,l,'Ground',[0   0 -1 min(t.Vertices(:,2))]);
s = add_shadow(t,l,'Ground',[0   0 -1 min(V(:,2))]);
D = matrixnormalize(s.Vertices(:,1)-s.Vertices(:,2));
s.FaceVertexCData = bsxfun(@times,D*0.3+0.7,[1 1 1]);
s.FaceColor = 'interp';
apply_ambient_occlusion(t,'AddLights',false);
%t.FaceVertexCData = bsxfun(@plus,0.5*t.FaceVertexCData,0.5*gray);
set(gca,'Visible','off')
axis equal;
view(4,16);
camproj('persp');
axis equal;
view(4,16);
camproj('persp');
set(gcf,'Color','w');

imwrite(myaa('raw'),'squiggly-torus.png');
%imwrite(myaa('raw'),'csaszar-torus.png');
