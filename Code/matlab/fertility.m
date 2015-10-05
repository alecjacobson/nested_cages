!rm .shrink_fine_expand_coarse_3D.cache.*.mat
load('~/Downloads/fert_initial_decimations (1).mat');
[V0,F0] = load_mesh('../../Meshes/Results/fert_varap_final/fert_0.obj');
levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
%V0 = V0*R;
%for c = numel(V_coarse)
%  V_coarse{c} = V_coarse{c}*R;
%end

render_data = RenderData();
render_data.R = axisangle2matrix([1 0 0],-pi/2);
render_data.prefix = 'fertility';
render_data.view = [-162,16];
render_data.ground = -0.49253;
render_data.axis = [-0.3 1.2 -0.3 1 render_data.ground-1e-2 0.7];
render_data.update_ao = true;
render_data.update_count = 0;
render_data.layer = 1;
render_data.zoom = 1.4;
render_data.colors = [
  244 243 73
  112 205 214
  213 122 213
  209 60 61
  98 180 71 
  36 39 201
  143 178 74]/255;
render_data.aa = false;

clf;
hold on;
render_data.tc = trisurf([1 1 1],0,0,0);
render_data.tf = trisurf([1 1 1],0,0,0,'EdgeColor','none');
render_data.l = light('Position',[-0.2 0.2 1],'Style','infinite');
camproj('persp');
[render_data.sc,~,render_data.csm] = add_shadow( ...
  render_data.tc, ...
  render_data.l,'Ground',[0 0 -1 render_data.ground-5e-3]);
render_data.sc.FaceColor = [1 1 1]*0.9;
[render_data.sf,~,render_data.fsm] = add_shadow( ...
  render_data.tf, ...
  render_data.l,'Ground',[0 0 -1 render_data.ground]);
render_data.sf.FaceColor = [1 1 1]*0.8;
hold off;
view(render_data.view(1),render_data.view(2));
axis equal;
axis manual;
axis(render_data.axis);
camzoom(render_data.zoom);
set(gca,'Visible','off')
set(gcf,'color','w');
set(0,'defaultaxesposition',[0 0 1 1]);
T = get(gca,'tightinset');
set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]);


[cages_V,cages_F,~,~,~,timing] = ...
multires_per_layer( ...
V0,F0, ...
true(1,numel(V_coarse)),'V_coarse',V_coarse,'F_coarse',F_coarse, ...
'QuadratureOrder',2, ...
'ExpansionEnergy','volumetric_arap', ...
'FinalEnergy','none', ...
'BetaInit',1e-2, ...
'EpsFinal',1e-3,...
'RenderData',render_data, ...
'EpsExpansion',1e-3);
