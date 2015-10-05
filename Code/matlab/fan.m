[OV,F] = load_mesh('fan.obj');
h = 1;
OV = [OV;[0 0 -h]];
O = outline(F);
F = [F;ones(size(O,1),1)*size(V,1) fliplr(O)];

CM = [
141,211,199
255,255,179
190,186,218
251,128,114
128,177,211
253,180,98
179,222,105
252,205,229
217,217,217
188,128,189
204,235,197
255,237,111]/255;
C = mod(1:size(F,1),12)';
colormap(CM);

OV(OV(:,3)==0.5,3) = h;
OV(OV(:,3)==h,1:2) = OV(OV(:,3)==h,1:2)*3.8
t = tsurf(F,OV,'FaceAlpha',0.8,'LineWidth',3,'CData',C);

l = light('Position',[-0.8 -0.9 1.0],'Style','infinite');
apply_ambient_occlusion(t,'AddLights',false);
camproj('persp');
axis equal;
set(gca,'Visible','off')
set(gcf,'color','w');
view(-115,36);


for pass = 1:2
  
  s = add_shadow(t,l,'Ground',[0   0 -1 min(OV(:,3))]);
  D = matrixnormalize(sum(bsxfun(@times,s.Vertices(:,1:2),l.Position(1:2)),2));
  s.FaceVertexCData = bsxfun(@times,D*0.1+0.9,[1 1 1]);
  s.FaceColor = 'interp';
  
  if pass == 1
    axis manual;
  end
   
  imwrite(myaa('raw'),sprintf('vouga-fan-%d.png',pass));

  delete(s);
  OV(OV(:,3)==0,:) = OV(OV(:,3)==0,:)*axisangle2matrix([0 0 1],-pi/3);
  t.Vertices = OV;
  t.CData = C;
end

%
%V1 = V;
%V1(OV(:,3)==0.0,1:2) = V1(OV(:,3)==0.0,1:2)*3.47
%V1(OV(:,3)==0.0,:) = V1(OV(:,3)==0.0,:)*axisangle2matrix([0 0 1],-pi/3);
%V1(OV(:,3)==0.0,3) = h;
%tsurf(F,V1,'FaceAlpha',0.6,'LineWidth',3,'CData',C);
%axis equal;
%pause
%
%V2 = V;
%V2(OV(:,3)==0.0,1:2) = V2(OV(:,3)==0.0,1:2)*3.47;
%tsurf(F,V2,'FaceAlpha',0.6,'LineWidth',3,'CData',C);
%axis equal;
%pause
%
%V3 = OV;
%V3(OV(:,3)==h,1:2) = V3(OV(:,3)==h,1:2)*3.47;
%tsurf(F,V3,'FaceAlpha',0.6,'LineWidth',3,'CData',C);
%axis equal;
