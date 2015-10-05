
[V0,F] = load_mesh('../../Meshes/Results/sca09_comparisons/horse_gallop_arap_0.obj');
[CV,CF] = load_mesh('../../Meshes/Results/sca09_comparisons/horse_gallop_arap_1.obj');
[CVs,CFs] = load_mesh('../../Meshes/Results/sca09_comparisons/horse-gallop-reference_cage.off');
[CVp,CFp] = load_mesh('../../Meshes/Results/ProgHulls/gallop_arap/PH_single/horse_001_1.obj');

R = axisangle2matrix([1 0 0],-pi/2);
CV = CV*R;
CVs = CVs*R;
CVp = CVp*R;
V0 = V0*R;

clf;
data = [];
T = get(gca,'tightinset');set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
for pass = 1:4 
  BB = [min([CV;CVs;CVp]);max([CV;CVs;CVp])];
  switch pass
  case 1
    data = render_in_cage(V0,F,CVs,CFs,'BoundingBox',BB,'View',[-66,18],'Data',data,'LightPosition',[-0.8 -0.9 0.8],'AmbientOcclusion',true);
  case 2
    data = render_in_cage(V0,F,CV,CF  ,'BoundingBox',BB,'View',[-66,18],'Data',data,'LightPosition',[-0.8 -0.9 0.8],'AmbientOcclusion',true);
    drawnow;
    drawnow;
  case 3
    data = render_in_cage(V0,F,CVp,CFp,'BoundingBox',BB,'View',[-66,18],'Data',data,'LightPosition',[-0.8 -0.9 0.8],'AmbientOcclusion',true);
  case 4
    delete(data.tc);
    delete(data.sc);
    delete(data.cc);
  end
  data.t.FaceAlpha = 1.0;
  %drawnow
  imwrite(myaa('raw'),sprintf('horse-gallop-%d.png',pass));
end
