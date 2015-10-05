%[V0,F0] = load_mesh('../../Meshes/Results/octopus-500k-volume_final/octopus_0.obj');
level = 2000;
[cages_V,cages_F,Pall] = multires_per_layer(V0,F0,level,'QuadratureOrder',2,'ExpansionEnergy','volumetric_arap','FinalEnergy','none','Smoothing',1e-1);
%[CV,CF] = decimate_cgal(V0,F0,level/size(F0,1));
%[CV,CF] = readPLY('octopus-2000.ply');
[V2] = load_mesh('octopus-rigid.obj');
[R,t] = fit_rigid(V2,V0);
clf;


cen0 = centroid(V0,F0);
cen2 = centroid(V2,F0);

CV2 = bsxfun(@plus,CV*R,t_new);
num_ip = size(intersect_other(CV,CF,CV2,CF),1)

dirs = cen0-cen2;
%dirs(2,:) = normalizerow(cross([1 0 0],dirs(1,:)))*norm(dirs(1,:));
%dirs(3,:) = normalizerow(cross(dirs(1,:),dirs(2,:)))*norm(dirs(1,:));

alpha = 1;

beta = 1;
while beta > 1e-6
  t_new = t+beta*(cen0-cen2);
  R_new = R*axisangle2matrix([1 0 0],beta*pi);
  CV2_new = bsxfun(@plus,CV*R,t_new);
  num_ip = size(intersect_other(CV,CF,CV2_new,CF),1);
  if num_ip > 0
    beta = beta/2;
    continue;
  end
  t = t_new;
  CV2 = bsxfun(@plus,CV*R,t);
  
  clf;
  hold on;
  %t1 =   tsurf(F0,V0,'FaceAlpha',alpha,'EdgeAlpha',alpha,'FaceColor','r','EdgeColor','none');
  %t2 =   tsurf(F0,V2,'FaceAlpha',alpha,'EdgeAlpha',alpha,'FaceColor','b','EdgeColor','none');
  tc1 =  tsurf(CF,CV,'FaceAlpha',alpha,'EdgeAlpha',alpha,'FaceColor','r');
  tc2 = tsurf(CF,CV2,'FaceAlpha',alpha,'EdgeAlpha',alpha,'FaceColor','b');
  camlight;
  hold off;
  view(95,10);
  axis equal;
  drawnow;
end

P = intersect_other(V0,F0,V2,F0);
F0I = F0(unique(P(:,1)),:);
conncomp(F0I)
