% Examples that worked with previous code (checking after re-writing of the
% code)
% 1) Flow step: s = 1e-3;
% 2) beta_init = 1e-2;
% 3) eps_proximity = 1e-4;

% [V0,F0] = readOBJ('../../Meshes/Results/anchor_volume/anchor.obj');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'eps_distance',1e-3);
% save('../../Meshes/Results/gargo_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK, using eps = 1e-3

[V0,F0] = readOBJ('../../Meshes/Results/couplingdown_volume/couplingdown.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'eps_distance',1e-3);
save('../../Meshes/Results/couplingdown_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.:

[V0,F0] = readOBJ('../../Meshes/Results/gargo_volume/gargo.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))]);
save('../../Meshes/Results/gargo_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.:

[V0,F0] = readOBJ('../../Meshes/Results/pelvis_volume/pelvis.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))]);
save('../../Meshes/Results/pelvis_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.:

[V0,F0] = readOBJ('../../Meshes/Results/Model3_volume/Model3.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))]);
save('../../Meshes/Results/Model3_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.:
