% % % Model9 - eps = 1e-3
% [V0,F0] = load_mesh('Model9_0.obj');
% [V_coarse{1},F_coarse{1}] = load_mesh('Model9_input_1.obj');
% [V_coarse{1},F_coarse{1}] = meshfix(V_coarse{1},F_coarse{1});
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, [0],...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_path', ...
% 'FinalEnergy','volume', ...
% 'BetaInit',1e-2, ...
% 'EpsExpansion',1e-2,...
% 'EpsFinal',1e-3,...
% 'V_coarse',V_coarse,'F_coarse',F_coarse,...
% 'PartialPath','partial_01_16.mat');
% write_cages('Model9_0001',cages_V,cages_F);
% 
% % Model9 - eps = 2e-4
% [V0,F0] = load_mesh('Model9_0.obj');
% [V_coarse{1},F_coarse{1}] = load_mesh('Model9_input_1.obj');
% [V_coarse{1},F_coarse{1}] = meshfix(V_coarse{1},F_coarse{1});
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, [0],...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_path', ...
% 'FinalEnergy','volume', ...
% 'BetaInit',1e-2, ...
% 'EpsExpansion',2e-3,...
% 'EpsFinal',2e-4,...
% 'V_coarse',V_coarse,'F_coarse',F_coarse,...
% 'PartialPath','partial_01_16.mat');
% write_cages('Model9_00002',cages_V,cages_F);
% 
% % Model9 Multires - eps = 1e-3
% [V0,F0] = load_mesh('multires/Model9_0.obj');
% [V_coarse{1},F_coarse{1}] = load_mesh('multires/Model9_input_7.obj');
% [V_coarse{1},F_coarse{1}] = meshfix(V_coarse{1},F_coarse{1});
% [V_coarse{2},F_coarse{2}] = load_mesh('multires/Model9_input_6.obj');
% [V_coarse{2},F_coarse{2}] = meshfix(V_coarse{2},F_coarse{2});
% [V_coarse{3},F_coarse{3}] = load_mesh('multires/Model9_input_5.obj');
% [V_coarse{3},F_coarse{3}] = meshfix(V_coarse{3},F_coarse{3});
% [V_coarse{4},F_coarse{4}] = load_mesh('multires/Model9_input_4.obj');
% [V_coarse{4},F_coarse{4}] = meshfix(V_coarse{4},F_coarse{4});
% [V_coarse{5},F_coarse{5}] = load_mesh('multires/Model9_input_3.obj');
% [V_coarse{5},F_coarse{5}] = meshfix(V_coarse{5},F_coarse{5});
% [V_coarse{6},F_coarse{6}] = load_mesh('multires/Model9_input_2.obj');
% [V_coarse{6},F_coarse{6}] = meshfix(V_coarse{6},F_coarse{6});
% [V_coarse{7},F_coarse{7}] = load_mesh('multires/Model9_input_1.obj');
% [V_coarse{7},F_coarse{7}] = meshfix(V_coarse{7},F_coarse{7});
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, [0 0 0 0 0 0 0],...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_path', ...
% 'FinalEnergy','volume', ...
% 'BetaInit',1e-2, ...
% 'EpsExpansion',1e-2,...
% 'EpsFinal',1e-3,...
% 'V_coarse',V_coarse,'F_coarse',F_coarse,...
% 'PartialPath','partial_01_16.mat');
% write_cages('multires/Model9_0001',cages_V,cages_F);

% Model9 Multires - eps = 2e-4
[V0,F0] = load_mesh('multires/Model9_0.obj');
[V_coarse{1},F_coarse{1}] = load_mesh('multires/Model9_input_7.obj');
[V_coarse{1},F_coarse{1}] = meshfix(V_coarse{1},F_coarse{1});
[V_coarse{2},F_coarse{2}] = load_mesh('multires/Model9_input_6.obj');
[V_coarse{2},F_coarse{2}] = meshfix(V_coarse{2},F_coarse{2});
[V_coarse{3},F_coarse{3}] = load_mesh('multires/Model9_input_5.obj');
[V_coarse{3},F_coarse{3}] = meshfix(V_coarse{3},F_coarse{3});
[V_coarse{4},F_coarse{4}] = load_mesh('multires/Model9_input_4.obj');
[V_coarse{4},F_coarse{4}] = meshfix(V_coarse{4},F_coarse{4});
[V_coarse{5},F_coarse{5}] = load_mesh('multires/Model9_input_3.obj');
[V_coarse{5},F_coarse{5}] = meshfix(V_coarse{5},F_coarse{5});
[V_coarse{6},F_coarse{6}] = load_mesh('multires/Model9_input_2.obj');
[V_coarse{6},F_coarse{6}] = meshfix(V_coarse{6},F_coarse{6});
[V_coarse{7},F_coarse{7}] = load_mesh('multires/Model9_input_1.obj');
[V_coarse{7},F_coarse{7}] = meshfix(V_coarse{7},F_coarse{7});
[cages_V,cages_F,~,~,~,timing] = ...
multires_per_layer( ...
V0,F0, [0 0 0 0 0 0 0],...
'QuadratureOrder',2, ...
'ExpansionEnergy','displacement_path', ...
'FinalEnergy','volume', ...
'BetaInit',1e-2, ...
'EpsExpansion',2e-3,...
'EpsFinal',2e-4,...
'V_coarse',V_coarse,'F_coarse',F_coarse,...
'PartialPath','partial_01_16.mat');
write_cages('multires/Model9_00002',cages_V,cages_F);
