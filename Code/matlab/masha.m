[V,F] = load_mesh('masha/masha-body.ply');
[V,IM] = remove_unreferenced(V,F);
F = IM(F);
[U,G] = load_mesh('masha/masha-clothes.ply');
[U,IM] = remove_unreferenced(U,G);
G = IM(G);
F = [F;size(V,1)+G];
V = [V;U];
C = connected_components(F);
I = randperm(max(C));
tsurf(F,V,'CData',I(C),'EdgeColor','none',fphong)
axis equal;
CM = parula(15);
colormap(CM(1:9,:));
%                               num_faces: 36135
%                            num_vertices: 19117
%                               num_edges: 55156
%                num_connected_components: 110
%                      num_boundary_loops: 104
%                     num_small_triangles: 2570
%                        num_small_angles: 103
%                      num_close_vertices: 25
%              num_selfintersecting_pairs: 12872

[W,BC,DV,Q] = voxelize(V,F,320,'Pad',2,'Interior',false);
save('masha-voxel-320.mat','W','BC','DV','Q','V','F');
% Get triangles from quads
DF = [Q(:,[1 2 3]);Q(:,[1 3 4])];
% Remove unreferenced corners
[SV,IM] = remove_unreferenced(DV,DF);
% re-index
SF = IM(DF);
SQ = IM(Q);
C = connected_components(SF);
trisurf(SQ,SV(:,1),SV(:,2),SV(:,3),'CData',C);
SQ1 = SQ(C(SQ(:,1))==1,:);
[SV1,IM] = remove_unreferenced(SV,SQ1);
SQ1 = IM(SQ1);
SF1 = SF(C(SF(:,1))==1,:);
SF1 = IM(SF1);
writePLY('masha-voxel-quads.ply',SV1,SQ1);
writePLY('masha-voxel-tris.ply',SV1,SF1);
[mV,mF] = meshfix(SV,SF1);

W = imfill(W,'holes');
FV = isosurface(W,0.5);
[mV,mF] = meshfix(FV.vertices,FV.faces);
writePLY('masha-voxel-contour-meshfix.ply',mV,mF);

W = imfill(W,'holes');
FV = isosurface(W,0.5);tsurf(FV.faces,FV.vertices);axis equal;
[mV,mF] = meshfix(FV.vertices,FV.faces);
writePLY('masha-voxel-contour-meshfix.ply',mV,mF);
writeOBJ('masha-voxel-contour-meshfix.obj',mV,mF);
[V,F] = load_mesh('masha-voxel-contour-meshfix-open-flipper.obj');
[mV,mF] = meshfix(V,F);

[mV,mF] = load_mesh('masha-voxel-contour-meshfix.ply');
mV = bsxfun(@plus,bsxfun(@times,bsxfun(@rdivide,mV-1, ...
  [size(W,2) size(W,1) size(W,3)]-1),max(BC)-min(BC)),min(BC));
[CV,CF] = decimate_cgal(mV,mF,36000/size(mF,1)*1.5);
[cages_V,cages_F,Pall] = multires_per_layer(V,F,[true],'V_coarse',{CV},'F_coarse',{CF},'QuadratureOrder',2,'ExpansionEnergy','volumetric_arap','FinalEnergy','none','BetaInit',1,'Eps',1e-4,'ExpandEvery',10);

