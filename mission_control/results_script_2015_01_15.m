% [V0,F0] = load_mesh('../../Meshes/Results/anchor_volume/anchor_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','surface_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'Eps',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/anchor_varap_new/anchor',cages_V,cages_F);
% save('../../Meshes/Results/anchor_varap_new/timing.mat','timing')
% % Obs.: VARAP --> Simulation stuck in 1st layer/last step
% % ARAP --> Taking too many iterations

% [V0,F0] = load_mesh('../../Meshes/Results/anchor_volume/anchor_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','displacement_path', ...
%   'FinalEnergy','volume', ...
%   'BetaInit',1e-2, ...
%   'Eps',5e-4,...
%   'EpsExpansion',5e-3,...
%   'PartialPath','partial_01_15.mat',...
%   'NewFlow',true);
% write_cages('../../Meshes/Results/anchor_volume_new_test/anchor',cages_V,cages_F);
% save('../../Meshes/Results/anchor_volume_new_test/timing.mat','timing')
% % Obs.: OK, 5 layers

% [V0,F0] = load_mesh('../../Meshes/Results/Model3_volume/Model3_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','displacement_path', ...
%   'FinalEnergy','volume', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',5e-4,...
%   'EpsExpansion',5e-3,...
%   'PartialPath','partial_01_15.mat',...
%   'NewFlow',true);
% write_cages('../../Meshes/Results/Model3_volume_new_test/Model3',cages_V,cages_F);
% save('../../Meshes/Results/Model3_volume_new_test/timing.mat','timing')

% [V0,F0] = load_mesh('../../Meshes/Results/arma_volume/arma_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat',...
%   'NewFlow',false);
% write_cages('../../Meshes/Results/arma_varap_new_no-opt-bug/arma',cages_V,cages_F);
% save('../../Meshes/Results/arma_varap_new_no-opt-bug/timing.mat','timing')

% [V0,F0] = load_mesh('../../Meshes/Results/arma_volume/arma_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-5, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat',...
%   'NewFlow',false,...
%   'BETA_MIN',1e-7);
% write_cages('../../Meshes/Results/arma_varap_new_no-opt-bug_beta-small/arma',cages_V,cages_F);
% save('../../Meshes/Results/arma_varap_new_no-opt-bug_beta-small/timing.mat','timing')
% % Obs.: It looked as bad as before (decreasing beta didn't help)

% [V0,F0] = load_mesh('../../Meshes/Results/Model3_volume/Model3_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/Model3_varap_final/Model3',cages_V,cages_F);
% save('../../Meshes/Results/Model3_varap_final/timing.mat','timing')
% % Obs.: 
% 
% [V0,F0] = load_mesh('../../Meshes/Results/Model7_volume/Model7_0.obj');
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','displacement_path', ...
%   'FinalEnergy','volume', ...
%   'BetaInit',5e-3, ...
%   'EpsFinal',5e-4,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/Model7_volume_final/Model7',cages_V,cages_F);
% save('../../Meshes/Results/Model7_volume_final/timing.mat','timing')
% % Obs.: 

% % New Alec's meshes
% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo-2/alien.ply');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/alien_varap_final/alien',cages_V,cages_F);
% save('../../Meshes/Results/alien_varap_final/timing.mat','timing')
% % Obs.: Simulation got stuck in the last steps of the sixth layer

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo-2/animal.obj');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','surface_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/animal_arap_final/animal',cages_V,cages_F);
% save('../../Meshes/Results/animal_arap_final/timing.mat','timing')
% % Obs.: Varap: Tetgen readFACE error (probably some internal error).
% % ARAP: OK!

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo-2/igea.obj');
% V0 = V0/max(max(abs(V0)));
% F0 = [F0(:,1) F0(:,3) F0(:,2)]; % inverted triangles
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/igea_varap_final/animal',cages_V,cages_F);
% save('../../Meshes/Results/igea_varap_final/timing.mat','timing')
% % Obs.: Input mesh is bad

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo-2/lobster-31k.ply');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',3, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/lobster_varap_final/animal',cages_V,cages_F);
% save('../../Meshes/Results/lobster_varap_final/timing.mat','timing')
% % Obs.: OK, 5 layers

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bimba-R8500.off');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',3, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/bimba_varap_final/bimba',cages_V,cages_F);
% save('../../Meshes/Results/bimba_varap_final/timing.mat','timing')
% % Obs.: OK

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bimba-R8500.off');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',3, ...
%   'ExpansionEnergy','surface_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/bimba_arap_final/bimba',cages_V,cages_F);
% save('../../Meshes/Results/bimba_arap_final/timing.mat','timing')
% % Obs.:
% 
% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bimba-R8500.off');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',3, ...
%   'ExpansionEnergy','displacement_path', ...
%   'FinalEnergy','volume', ...
%   'BetaInit',1e-2, ...
%   'EpsExpansion',5e-3,...
%   'EpsFinal',5e-4,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/bimba_volume_final/bimba',cages_V,cages_F);
% save('../../Meshes/Results/bimba_volume_final/timing.mat','timing')
% % Obs.:
% 
% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bimba-R8500.off');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',3, ...
%   'ExpansionEnergy','displacement_initial', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsExpansion',1e-3,...
%   'EpsFinal',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/bimba_dispinitial_final/bimba',cages_V,cages_F);
% save('../../Meshes/Results/bimba_dispinitial_final/timing.mat','timing')
% % Obs.:

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bunny-50k.ply');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/bunny-50k_varap_final/bunny',cages_V,cages_F);
% save('../../Meshes/Results/bunny-50k_varap_final/timing.mat','timing')
% % Obs.: OK

% % Re-run to try to obtain 5 layers
% [V0,F0] = load_mesh('../../Meshes/Results/buhda-100k-clean/buhda-100k.ply');
% % Had to clean up to help the flow
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
%   multires_per_layer( ...
%   V0,F0, ...
%   levels, ...
%   'QuadratureOrder',2, ...
%   'ExpansionEnergy','volumetric_arap', ...
%   'FinalEnergy','none', ...
%   'BetaInit',1e-2, ...
%   'EpsFinal',1e-3,...
%   'EpsExpansion',1e-3,...
%   'ExpandEvery', 10,...
%   'PartialPath','partial_01_15.mat');
% write_cages('../../Meshes/Results/buhda_varap_inflated/buhda',cages_V1,cages_F1);
% save('../../Meshes/Results/buhda_varap_inflated/timing.mat','timing')
% % Obs.: 5 layers, but inflated

% [V0,F0] = load_mesh('../../Meshes/Results/buhda-100k-clean/buhda-100k.ply');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, ...
% levels, ...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_initial', ...
% 'FinalEnergy','none', ...
% 'BetaInit',5e-1, ...
% 'EpsFinal',1e-3,...
% 'EpsExpansion',1e-3,...
% 'ExpandEvery', 2,...
% 'PartialPath','partial_01_15.mat',...
% 'BETA_MIN',5e-3);
% write_cages('../../Meshes/Results/buhda_dispinitial_beta05_expand2/buhda',cages_V,cages_F);
% save('../../Meshes/Results/buhda_dispinitial_beta05_expand2/timing.mat','timing')
% % Obs.: Stuck at the end of the second layer

% [V0,F0] = load_mesh('../../Meshes/Results/buhda-100k-clean/buhda-100k.ply');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, ...
% levels, ...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_initial', ...
% 'FinalEnergy','none', ...
% 'BetaInit',5e-2, ...
% 'EpsFinal',1e-3,...
% 'EpsExpansion',1e-3,...
% 'ExpandEvery', 2,...
% 'PartialPath','partial_01_15.mat',...
% 'BETA_MIN',5e-4);
% write_cages('../../Meshes/Results/buhda_dispinitial_beta005_expand2/buhda',cages_V,cages_F);
% save('../../Meshes/Results/buhda_dispinitial_beta005_expand2/timing.mat','timing')
% % Obs.: Stuck at the end of the second layer

% [V0,F0] = load_mesh('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model5_fixed.off');
% V0 = V0/max(max(abs(V0)));
% levels = floor(2.^((-14:2:-2)/3)*size(F0,1));
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, ...
% levels, ...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_initial', ...
% 'FinalEnergy','none', ...
% 'BetaInit',1e-2, ...
% 'EpsFinal',1e-3,...
% 'EpsExpansion',1e-3,...
% 'ExpandEvery', 2,...
% 'PartialPath','partial_01_15.mat',...
% 'BETA_MIN',1e-4);
% write_cages('../../Meshes/Results/Model5_dispinitial_beta001_expand2/Model5',cages_V,cages_F);
% save('../../Meshes/Results/Model5_dispinitial_beta001_expand2/timing.mat','timing')
% % Obs.: Generated 1 layer. It would take too long to minimize 'disp_initial' per step

% % mv partial_01_15.mat partial_01_15_Model5_1.mat 
% load('partial_01_15_Model5_1.mat');
% V0 = cages_V{7};
% F0 = cages_F{7};
% levels = floor(2.^((-10:2:-2)/3)*size(F0,1)); %should be one more layer
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, ...
% levels, ...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_path', ...
% 'FinalEnergy','displacement_initial', ...
% 'BetaInit',1e-2, ...
% 'EpsFinal',1e-3,...
% 'EpsExpansion',1e-3,...
% 'ExpandEvery', 2,...
% 'PartialPath','partial_01_15.mat',...
% 'BETA_MIN',1e-4);
% write_cages('../../Meshes/Results/Model5_dispinitial_beta001_expand2/Model5',cages_V,cages_F);
% save('../../Meshes/Results/Model5_dispinitial_beta001_expand2/timing.mat','timing')
% % Generated 1 layer and crashed because coarse mesh self-inetrsects

% % mv partial_01_15.mat partial_01_15_Model5_2.mat 
% load('partial_01_15_Model5_2.mat');
% V0 = cages_V{5};
% F0 = cages_F{5};
% levels = floor(2.^((-10:2:-2)/3)*size(F0,1))
% levels(5) = levels(5)+1000; %trick to avoid self-intersections
% [cages_V,cages_F,~,~,~,timing] = ...
% multires_per_layer( ...
% V0,F0, ...
% levels, ...
% 'QuadratureOrder',2, ...
% 'ExpansionEnergy','displacement_path', ...
% 'FinalEnergy','displacement_initial', ...
% 'BetaInit',1e-2, ...
% 'EpsFinal',1e-3,...
% 'EpsExpansion',1e-3,...
% 'ExpandEvery', 2,...
% 'PartialPath','partial_01_15.mat',...
% 'BETA_MIN',1e-4);
% write_cages('../../Meshes/Results/Model5_dispinitial_beta001_expand2/Model5',cages_V,cages_F);
% save('../../Meshes/Results/Model5_dispinitial_beta001_expand2/timing.mat','timing')
% % Stuck in simulation last step of second level (at least got one more)

% mv partial_01_15.mat partial_01_15_Model5_3.mat 
load('partial_01_15_Model5_3.mat');
V0 = cages_V{5};
F0 = cages_F{5};
levels = floor(2.^((-8:2:-2)/3)*size(F0,1))
[cages_V,cages_F,~,~,~,timing] = ...
multires_per_layer( ...
V0,F0, ...
levels, ...
'QuadratureOrder',2, ...
'ExpansionEnergy','displacement_path', ...
'FinalEnergy','displacement_initial', ...
'BetaInit',1e-2, ...
'EpsFinal',1e-3,...
'EpsExpansion',1e-3,...
'PartialPath','partial_01_15.mat',...
'BETA_MIN',1e-4);
write_cages('../../Meshes/Results/Model5_dispinitial_beta001_expand2/Model5',cages_V,cages_F);
save('../../Meshes/Results/Model5_dispinitial_beta001_expand2/timing.mat','timing')

% Other models for the zoo
% 1) handles smooth
% 2) casting-10k (refine to 40k)
% 3) rocker arm
% 4) SGP/davidone-s20k
% 5) SGP/chinese-dragon
% 6) SGP/ramesses-20k