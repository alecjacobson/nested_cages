% results generated with:
% 1) Flow step: s = 1e-3;
% 2) beta_init = 1e-1;
% 3) eps_proximity = 5e-4;
% 4) Eltopo for inflation with 10*eps_proximity

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/anchor-11k.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/anchor_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: Last level lost too many details

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/couplingdown_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/couplingdown_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: Last levels lost too many details
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/pelvis-40k-butterfly.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/pelvis_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/gargo-R13k.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/gargo_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/horse_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/horse_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK
% 
% [V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/l_brachioradialis_2k.obj');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-5,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/l_brachioradialis_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK. Used eps=5e-5 to get a less bumpy result. For coarse layers,
% % flow has bumpiness and this is passed to the cage via inflation if eps
% % is high.

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model1_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model1_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model3_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model3_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model7_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model7_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model9_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model9_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model4_fixed.off');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model4_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.: Letting this for the end, harder simulation (knees too close)

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model5_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model5_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: Generated 2 layers, but difficult flow (let for later)
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model6_fixed.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/Model6_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: Difficult flow, let for later

% Saving results
% cd '/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/'
% load('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 -
% multires_results/pelvis_volume.mat');
% mkdir pelvis_volume
% writeOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 -
% multires_results/pelvis_volume/pelvis.obj',V0,F0);
% for k=1:size(cages_V,2)-1
% filename = sprintf('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/pelvis_volume/pelvis_%d.obj',size(cages_V,2)-k);
% writeOBJ(filename,cages_V{k},cages_F{k});
% end