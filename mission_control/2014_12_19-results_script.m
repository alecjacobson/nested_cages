% results generated with:
% 1) Flow step: s = 1e-3;
% 2) eps_proximity = 5e-4
% 3) Second order quadrature
% 4) Per layer volume minimization
% 5) No armijo seacrh for flow step

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/pelvis-4k.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/pelvis_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK. Generate with less gaps to have more layers.
% 
% [V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/l_brachioradialis_2k.obj');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/l_brachioradialis_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK. Generate with less gaps to have more layers.
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/anchor-11k.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/anchor_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK. Generate with less gaps to have more layers.

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bunny-35k.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/bunny_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/gargo-R13k.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/gargo_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/horse_fixed.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/horse_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK
% 
% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/arma.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/arma_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK

% [V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/octopus_100k.obj');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/octopus_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/bunny-35k.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/bunny_finer_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/rampant.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/8) floor(size(F0,1)/4)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/rampant_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: It generated only 2 layers. Run with smaller gaps.

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/couplingdown_fixed.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/couplingdown_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: Coarser layers don't resamble original shape. Run with smaller
%   gaps.

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model1_fixed.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/Model1_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: OK, generated with w_lap = 0.1 (didn't work with w_lap=0.0)

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model2_fixed.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/Model2_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: w_lap = 0.1. At the end of the second layer ElTopo didn't work.
% % Velocity filter took very long inside Matlab and crashed outside it (seg fault 11).

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model3_fixed.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/Model3_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: At the end of the fifth layer, velocity filter crashed inside
% % Matlab. 

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/rampant.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/rampant_finer_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: It generated only 2 layers. Run with smaller gaps.

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/casting-10k.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/casting_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: Generated 2 layers then decimation was too different from original.

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model13_fixed.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/Model13_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: w_lap = 0.0 didn't work. Trying w_lap = 0.1. Didn't work either.
% Trying with s=5e-4. Didn't work either.

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/Model4_fixed.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/Model4_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: 

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/meshes-for-leo/beast_fixed.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-14/3)*floor(size(F0,1)) 2^(-12/3)*floor(size(F0,1)) 2^(-10/3)*floor(size(F0,1)) 2^(-8/3)*floor(size(F0,1)) 2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*size(F0,1))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/beast_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: w_lap = 0.0 didn't work for the fisrt layers.

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/dragon-S10K.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*(size(F0,1)))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/dragon_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: Expanding at every 50 iterations worked for the first layer

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/dragon_recon/dragon_30k.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[2^(-6/3)*floor(size(F0,1)) floor(2^(-4/3)*size(F0,1)) floor(2^(-2/3)*(size(F0,1)))],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/dragon_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: Rank deficient matrix when trying to use smoothing

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/heart_fixed.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/heart_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: Input mesh self-intersects. Removing intersections leaves no gap
% between some parts of the maseh (impossible to generate a cage). Also,
% flow didn't work

[V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/happy_recon/happy_100k.off');
tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/happy_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% Obs.: Complicated flow

% [V0,F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/iwires-alien-19k.off');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/iwires-alien_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: 

% [V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/microscope-14k.obj');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/microscope_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: Difficult flow

% [V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/heart_palmonary_artery_46k.obj');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/palmonary_artery_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: CGAL simplification didn't work on this one

% [V0,F0] = readOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/models/heart_leftatrium_70k.obj');
% tic; [cages_V,cages_F,Pall,V_coarse,F_coarse] = multires_per_layer(V0,F0,[floor(size(F0,1)/128) floor(size(F0,1)/64) floor(size(F0,1)/32) floor(size(F0,1)/16) floor(size(F0,1)/8) floor(size(F0,1)/4) floor(size(F0,1)/2)],'method','shrink_fine_and_expand_coarse','quadrature_order',2); toc;
% save('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/leftatrium_volume.mat','cages_V','cages_F','Pall','V_coarse','F_coarse','V0','F0');
% % Obs.: CGAL simplification didn't work on this one


% Saving results
% cd '/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/'
% load('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 -
% multires_results/pelvis_volume.mat');
% mkdir pelvis_volume
% writeOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 -
% multires_results/pelvis_volume/pelvis.obj',V0,F0);
% for k=1:size(cages_V,2)-1
% filename = sprintf('/Users/Leo/PHD_Work/Cage_Generation_2013/results/2014.12.19 - multires_results/pelvis_volume/pelvis_%d.obj',size(cages_V,2)-k);
% writeOBJ(filename,cages_V{k},cages_F{k});
% end