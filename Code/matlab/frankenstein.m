addpath ~/Dropbox/multires/
[V,F] = load_mesh('~/Dropbox/models/frankenstein-orient_outward_ao-oriented.obj');
[V0,F0] = load_mesh('../../Meshes/Results/frankenstein_varap_final/frankenstein_0.obj');
V = bsxfun(@minus,V,min(V))/max(max(V)-min(V));
V = bsxfun(@plus,V*max(max(V0)-min(V0)),min(V0));
[CV,CF] = load_mesh('../../Meshes/Results/frankenstein_varap_final/frankenstein_1.obj');
R = axisangle2matrix([0 1 0],pi/2);

%[V,F] = upsample(V,F);

CV = CV*R;
V = V*R;
V0 = V0*R;
C = connected_components(F);
%I = randperm(max(C));
I = [24 25 2 20 15 10 16 21 13 7 22 9 12 14 1 6 5 19 11 23 8 3 17 4 18];
C = I(C);
%clf;
%hold on;
%t = tsurf(F,V,'CData',C,'EdgeColor','none');
%tsurf(CF,CV,'FaceAlpha',0.2,'EdgeAlpha',0.2,'LineWidth',4);
%hold off;
%axis equal;
%camproj('persp');
%%apply_ambient_occlusion();
%
%h = mean(mean(edge_lengths(CV,CF)));
%max_vol = h^3;
%[DV,DT,DF] = tetgen(CV,CF,'Flags',sprintf('-a%0.17f -q1.5Y',max_vol));
%
%P = prolongation(DV,DT,V,'Extrapolation','linear');
%
%M = massmatrix(DV,DT);
%stiff = 1.2e-1;
%M = M*(stiff^-2);
%
%
%ground = -1.1;
%g = 1e-4;
%h = 1;
%
%s = [];
%sc = [];
%data = [];
%DU = DV;
%DUm1 = DU;
%DUm1 = DUm1*axisangle2matrix([0 0 1],-pi/300);
%DUm1 = bsxfun(@plus,bsxfun(@minus,DUm1,[0 0 mean(DU(:,3))])*axisangle2matrix([0 1 0],pi/700),[0 0 mean(DU(:,3))]);
%%tc = tsurf(DF,DV,'FaceAlpha',0.2,'EdgeAlpha',0.2,'LineWidth',4);
%clf;
%hold on;
%t = tsurf(F,P*DV,'CData',C,'EdgeColor','none');
%%tc = tsurf(DF,DV,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.4,'LineWidth',1);
%hold off;
%l = light('Position',[-0.3 -0.1 1 ],'Style','infinite');
%%[sc,~,dshadow_mat] = add_shadow(tc,l,'Ground',[0 0 -1 ground-1e-4]);
%sc.FaceColor = [1 1 1]*0.8;
%[s,~,shadow_mat] = add_shadow(t,l,'Ground',[0 0 -1 ground]);
%axis equal;
%camproj('persp');
%set(gca,'Visible','off');
%set(gcf,'color','w');
%T = get(gca,'tightinset');
%set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);
%%view(-80,6);
%view(-133,26);
%axis([-1 1.5 -1.5 1.2 -1.3 0.75]);
%damp = 1e-3;

%for iter = 368:734
%  D = ground-DU(:,3);
%  D(D<0) = 0;
%  DU(D>0,3) = ground;
%  G = -g*ones(size(DU,1),1);
%  fext = [zeros(size(DU,1),2) M*(G+D)];
%  DU_prev = DU;
%  [DU,data] = arap( ...
%    DV,DT,[],[], ...
%    'Energy','elements', ...
%    'V0',DU, ...
%    'Vm1',DUm1, ...
%    'Data',data, ...
%    'TimeStep', h*stiff, ...
%    'MaxIter',500, ...
%    'Dynamic',fext);
%  DUm1 = DU_prev;
%  % really regretting pass positions instead of velocities
%  %
%  % v = u0 - u1
%  % vd = v*d
%  % v*d = (u0 - u1)*d
%  % v*d = u0*d - u1*d
%  % v*d = u0*(1+d-1) - u1*d
%  % v*d = u0 + u0*(d-1) - u1*d
%  % v*d = u0 - -u0*(d-1) - u1*d
%  % v*d = u0 - (u1*d-u0*(d-1))
%  % 
%  DUm1 = DUm1+damp*(DU-DUm1);
%
%  if mod(iter,2) == 1
%    %tc.Vertices = DU;
%    U = P*DU;
%    t.Vertices = U;
%    if mod(iter,4) == 1
%      t.FaceColor = 'interp';
%      t.CData = C;
%      apply_ambient_occlusion(t,'AddLights',false);
%    end
%    SU = [U ones(size(U,1),1)]*shadow_mat(:,:,2)';
%    SU = bsxfun(@rdivide,SU(:,1:3),SU(:,4));
%    s.Vertices = SU;
%    %SDU = [DU ones(size(DU,1),1)]*dshadow_mat(:,:,2)';
%    %SDU = bsxfun(@rdivide,SDU(:,1:3),SDU(:,4));
%    %sc.Vertices = SDU;
%
%    drawnow;
%    frame = getframe(gcf);
%    im = frame.cdata;
%    imwrite(im,sprintf('frankenstein-no-cage-%06d.jpg',floor(iter/2)));
%  end
%end

clf;
hold on;
[V,F] = upsample(V,F);
[V,F] = upsample(V,F);
C = connected_components(F);
%I = randperm(max(C));
I = [24 25 2 20 15 10 16 21 13 7 22 9 12 14 1 6 5 19 11 23 8 3 17 4 18];
C = I(C);
t = tsurf(F,V,'CData',C,'EdgeColor','none');
apply_ambient_occlusion(t,'AddLights',false);
tc = tsurf(CF,CV,'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.4,'LineWidth',1);
hold off;
l = light('Position',[-0.3 -0.1 1 ],'Style','infinite');
%[SV,SF] = signed_distance_isosurface(V0,F0,'SignedDistanceType','pseudonormal','Level',0.029,'DistanceBound',0.01854,'RadiusBound',0.01854);
[SV,SF] = readPLY('frankenstein-cgal.ply');
SV = SV*R;
ground = min([V(:,3);CV(:,3);SV(:,3)]);
%[sc,~,dshadow_mat] = add_shadow(tc,l,'Ground',[0 0 -1 ground-1e-3]);
%sc.FaceColor = [1 1 1]*0.8;
%[s,~,shadow_mat] = add_shadow(t,l,'Ground',[0 0 -1 ground]);
axis equal;
camproj('persp');
set(gca,'Visible','off');
set(gcf,'color','w');
%view(-80,6);
view(-100,4);
set(gca,'Visible','off');
set(gcf,'Color','w');
pause

set(tc,'Vertices',SV,'Faces',SF);
SSV = [SV ones(size(SV,1),1)]*dshadow_mat(:,:,2)';
SSV = bsxfun(@rdivide,SSV(:,1:3),SSV(:,4));
%set(sc,'Vertices',SSV,'Faces',SF);
pause


set(tc,'Vertices',V0,'Faces',F0);
SV0 = [V0 ones(size(V0,1),1)]*dshadow_mat(:,:,2)';
SV0 = bsxfun(@rdivide,SV0(:,1:3),SV0(:,4));
%set(sc,'Vertices',SV0,'Faces',F0);
pause
