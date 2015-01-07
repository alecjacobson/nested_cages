% results generated with:
% 1) Flow step: s = 1e-3;
% 2) beta_init = 1e-1;
% 3) eps_proximity = 5e-4;
% 4) Eltopo for inflation with 10*eps_proximity
% 

% sanity check
[V0,F0] = readOBJ('../../Meshes/Results/gargo_volume/gargo.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[floor(2^(-14/3)*size(F0,1)) floor(2^(-12/3)*size(F0,1)) floor(2^(-10/3)*size(F0,1)) floor(2^(-8/3)*size(F0,1)) floor(2^(-6/3)*size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
save('../../Meshes/Results/gargo_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.: OK

% mini stress test
[V0,F0] = readOBJ('../../Meshes/Results/octopus_volume/octopus.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[floor(2^(-14/3)*size(F0,1)) floor(2^(-12/3)*size(F0,1)) floor(2^(-10/3)*size(F0,1)) floor(2^(-8/3)*size(F0,1)) floor(2^(-6/3)*size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))]);
save('../../Meshes/Results/octopus_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.: OK