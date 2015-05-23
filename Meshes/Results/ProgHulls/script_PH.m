% basename = 'alien_varap/PH_multi/alien';
% num_levels = 5;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'anchor_volume/PH_multi/anchor'
% num_levels = 6;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'animal_arap/PH_multi/animal'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% % Obs.: This was skipped because slowed my computer down.

% basename = 'arma_volumetric_arap/PH_multi/arma'
% num_levels = 6;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'bimba_varap/PH_multi/bimba'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'bunny-50k_varap/PH_multi/bunny'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'couplingdown_volume/PH_multi/couplingdown'
% num_levels = 3;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% % [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'disney1_varap/PH_multi/disney1'
% num_levels = 5;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'hand_varap/PH_multi/hand'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'horse_varap_25/PH_multi/horse'
% num_levels = 25;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = 1000:1000:25000;
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'lobster_varap/PH_multi/lobster'
% num_levels = 5;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 

% basename = 'manhead_varap/PH_multi/manhead'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'maxplank_varap_50/PH_multi/maxplank'
% num_levels = 50;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = 750:250:13000;
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'Model1_varap/PH_multi/Model1'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'Model3_varap/PH_multi/Model3'
% num_levels = 4;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'Model4_varap/PH_multi/Model4'
% num_levels = 4;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

% basename = 'Model7_volume/PH_multi/Model7'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);
% 
% basename = 'Model9_varap/PH_multi/Model9'
% num_levels = 7;
% [V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
% [V0,F0] = meshfix(V0,F0);
% levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
% [hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
% write_cages(basename,hulls_V_mex,hulls_F_mex);

basename = 'pelvis_volume/PH_multi/pelvis'
num_levels = 7;
[V0,F0] = load_mesh(sprintf('%s_0.obj',basename));
[V0,F0] = meshfix(V0,F0);
levels = floor(2.^((-2*num_levels:2:-2)/3)*size(F0,1));
[hulls_V_mex,hulls_F_mex] = progressive_hulls_mex(V0,F0,levels,0.05);
[hulls_V_mex,hulls_F_mex] = progressive_hulls_clean(hulls_V_mex,hulls_F_mex);
write_cages(basename,hulls_V_mex,hulls_F_mex);