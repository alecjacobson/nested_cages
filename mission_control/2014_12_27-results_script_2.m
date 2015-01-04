% results generated with:
% 1) Flow step: s = 1e-3;
% 2) beta_init = 1e-1;
% 3) eps_proximity = 5e-4;
% 4) Eltopo for inflation with 10*eps_proximity

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bunny-35k.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/bunny_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: OK.

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/arma.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[floor(size(F0,1)/64) floor(size(F0,1)/8)],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/arma_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: Too bumpy result.
% 
[V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/octopus_100k.obj');
[cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/octopus_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% Obs.: Big example. Let it for later. It was taking more than 8 hours for
% the first layer. It was too bumpy, try with smaller eps

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/rampant-100k.off');
% [cages_V,cages_F,Pall,V_coarse,F_coarse,timing] = multires_per_layer(V0,F0,[floor(size(F0,1)/144) floor(size(F0,1)/12)],'method','shrink_fine_and_expand_coarse','quadrature_order',2,'eps_distance',5e-4,'beta_init',1e-1);
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.27 - multires_results/rampant_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0','timing');
% % Obs.: [floor(size(F0,1)/256) floor(size(F0,1)/16)] didn't work (bad flow)

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