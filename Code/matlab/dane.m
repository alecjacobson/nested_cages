%[V,F] = load_mesh('~/Downloads/meshes_cages2/dane-reference.off');
%F = F(doublearea(V,F)>eps,:);
%[V,~,IM] = remove_duplicate_vertices(V,0);
%F = IM(F);
%
%method = 2;
%switch method
%case 1
%  [SV,SF,~,~,IM] = selfintersect(V,F);
%  F = IM(SF);
%  [V,IM] = remove_unreferenced(SV,F);
%  F = IM(F);
%  [TV,TT] = cdt(V,F,'TetgenFlags','-T1e-10Y');
%  BC = barycenter(TV,TT);
%  W = winding_number(V,F,BC);
%  I = abs(W)>0.5;
%  CT = TT(I,:);
%  CF = boundary_faces(CT);
%  [CV,IM] = remove_unreferenced(TV,CF);
%  CF = IM(CF);
%  G = outer_hull(CV,G);
%  [V,F] = meshfix(CV,G);
%case 2
%  [W,BC,DV,Q] = voxelize(V,F,40,'Boundary',true);
%  FV = isosurface(W,0.5);
%  CVf = bsxfun(@plus,bsxfun(@times,bsxfun(@rdivide,FV.vertices-1, ...
%    [size(W,2) size(W,1) size(W,3)]-1),max(BC)-min(BC)),min(BC));
%  CFf = FV.faces;
%  C = connected_components(CFf);
%  CFf = CFf(C(CFf(:,1))==1,:);
%  writeOBJ('dane-voxel-cage.obj',CVf,CFf);
%  pause('decimate in openflipper to 674 triangles and save into this obj');
%  [CV0,CF0] = load_mesh('dane-voxel-cage.obj');
%  [mV,mF] = meshfix(CV0,CF0);
%  V0 = bsxfun(@minus,V,min(V))/max(max(V)-min(V));
%  CV0 = bsxfun(@minus,mV,min(V))/max(max(V)-min(V));
%  CF = mF;
%end

addpath ../../Meshes/Results/sca09_comparisons/
[V0,F] = load_mesh('dane.ply');
F = fliplr(F);
[CV,CF] = load_mesh('dane-voxel-cage-opt-ref.ply');
[CVs,CFs] = load_mesh('dane-sca-cage.ply');
[CVp,CFp] = load_mesh('dane-fixed-with-tail-ph.ply');

R = axisangle2matrix([1 0 0],-pi/2);
CV = CV*R;
CVs = CVs*R;
CVp = CVp*R;
V0 = V0*R;

clf;
data = [];
T = get(gca,'tightinset');set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
for pass = [2 1 3 4]
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
  imwrite(myaa('raw'),sprintf('dane-%d.png',pass));
end
