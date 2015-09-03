%[V,F] =   readOBJ('../../Meshes/Results/KenshiHand_planarity/KenshiHandFilled_planarity_0.obj');
%[CV0,CQ] = readOBJ('../../Meshes/Results/KenshiHand_planarity/KenshiHandFilled_planarity_1_quad.obj','Quads',true);
%
%rot = axisangle2matrix([-0.2 1 0],pi/2) * axisangle2matrix([1 0 0],pi/4) * axisangle2matrix([0 1 0],-pi/10);
%V = V*rot;
%CV0 = CV0*rot;
%
%CV = planarize(CV0,CQ);
%tetgen_flags = '-q2a0.000001';
%[DQV,DQT,DQF] = encage(V,CV,CQ,'TetgenFlags',tetgen_flags);
%
%% Cage Edges
%CE = [CQ(:,[1 2]); CQ(:,[2 3]); CQ(:,[3 4]); CQ(:,[4 1])];
%CE = unique(sort(CE,2),'rows');
%DQW = hc(DQV,DQT,DQF,CV,CE);
%
%CDodd =  CQ(:,[1 2 3 1 3 4]);
%CDeven = CQ(:,[1 2 4 3 4 2]);
%R = rand(size(CQ,1),1)>0.5;
%CFodd =  reshape([ CDodd(R,:);CDeven(~R,:)]',3,size(CQ,1)*2)';
%CFeven = reshape([CDeven(R,:); CDodd(~R,:)]',3,size(CQ,1)*2)';
%[DoddV,DoddT,DoddF,b,bc] = encage(V,CV,CFodd,'TetgenFlags',tetgen_flags);
%DoddW = kharmonic(DoddV,DoddT,b,bc);
%
%[DevenV,DevenT,DevenF,b,bc] = encage(V,CV,CFeven,'TetgenFlags',tetgen_flags);
%DevenW = kharmonic(DevenV,DevenT,b,bc);
%
%sel = [14 16 20 40 68 84 101 47];

CU = load_mesh('../../Meshes/Results/KenshiHand_planarity/KenshiHand_posed.obj');
CU = CU*rot;

%for pass = 1:3
%
%  switch pass
%  case 1
%    VV = DQV;
%    FF = DQF;
%    CD = [1 1 1];
%    W = DQW;
%  case 2
%    VV = DoddV;
%    FF = DoddF;
%    CD = CFodd;
%    W = DoddW;
%  case 3
%    VV = DevenV;
%    FF = DevenF;
%    CD = CFeven;
%    W = DevenW;
%  end
%
%  clf;
%  hold on;
%  t = tsurf(FF,VV,'CData',sum(W(:,sel),2),fphong,'EdgeColor','none');
%  te = trisurf(CD,CV(:,1),CV(:,2),CV(:,3),'EdgeColor','k','FaceColor','none','LineWidth',2);
%  tq = trisurf(CQ,CV(:,1),CV(:,2),CV(:,3),'EdgeColor','k','FaceColor','none','LineWidth',5);
%  colormap(parula(256));
%  view(-6,12);
%  l = light('Position',[0.2 -0.2 1],'Style','infinite');
%  s = add_shadow(t,l,'Ground',[0 0 -1 min(VV(:,3))]);
%  apply_ambient_occlusion(t,'AddLights',false);
%  axis equal
%  hold off;
%  set(gca,'Visible','off');
%  set(gcf,'Color','w');
%  pause();
%
%
%  DU = W*CU;
%  U = DU(size(CV,1)+(1:size(V,1)),:);
%  clf;
%  hold on;
%  t = tsurf(F,U,'FaceColor',[0.0482,0.3651,0.8722],'EdgeColor','none'); 
%  te = trisurf(CD,CU(:,1),CU(:,2),CU(:,3),'EdgeColor','k','LineWidth',2,'FaceColor','none','FaceAlpha',0.3);
%  tq = trisurf(CQ,CU(:,1),CU(:,2),CU(:,3),'EdgeColor','k','FaceColor',[0.9 0.9 0.9],'LineWidth',5,'FaceAlpha',0.3);
%  view(-6,12);
%  l = light('Position',[0.2 -0.2 1],'Style','infinite');
%  ground = [0 0 -1 min(VV(:,3))];
%  s = add_shadow(t,l,'Ground',ground);
%  sc = add_shadow(tq,l,'Ground',ground-[0 0 0 1e-4]);
%  s.FaceColor = [1 1 1]*0.8;
%  sc.FaceColor = [1 1 1]*0.95;
%  apply_ambient_occlusion(t,'AddLights',false);
%  axis equal
%  hold off;
%  set(gca,'Visible','off');
%  set(gcf,'Color','w');
%  pause();
%end
%

for pass = 1:2
  switch pass
  case 1
    CU = load_mesh('/Users/ajx/Documents/nested_cages/Meshes/Results/KenshiHand_planarity/KenshiHandFilled_input.ply');
    CU = CU*rot;
    alpha = 0.9;
  case 2
    CU = CV;
    alpha = 1.0;
  end
  clf;
  hold on;
  t = tsurf(F,V,'FaceColor',[0.0482,0.3651,0.8722],'EdgeColor','none','FaceAlpha',alpha); 
  tc = trisurf(CQ,CU(:,1),CU(:,2),CU(:,3),'EdgeColor','k','FaceColor',[0.9 0.9 0.9],'LineWidth',5,'FaceAlpha',0.3);
  view(-6,12);
  l = light('Position',[0.2 -0.2 1],'Style','infinite');
  ground = [0 0 -1 min(VV(:,3))];
  s = add_shadow(t,l,'Ground',ground);
  sc = add_shadow(tc,l,'Ground',ground-[0 0 0 1e-4]);
  s.FaceColor = [1 1 1]*0.8;
  sc.FaceColor = [1 1 1]*0.95;
  apply_ambient_occlusion(t,'AddLights',false);
  axis equal
  hold off;
  set(gca,'Visible','off');
  set(gcf,'Color','w');
  pause();
end


%W = DW(size(CV,1)+(1:size(V,1)),:);
%writeTGF('kenshi-hand-cage.tgf',CV,CE);
%writeDMAT('kenshi-hand-cage.dmat',W);
%writeOBJ('kenshi-hand-tri.obj',V,F);
