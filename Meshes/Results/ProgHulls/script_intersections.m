% basename = 'fert_volume/PH_multi/fert'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'gallop_arap/PH_single/horse_001'
% num_levels = 1;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'gargo_volume/PH_multi/gargo'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'handles_volume/PH_single/handles'
% num_levels = 1;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'homer_volume/PH_multi/homer'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'homer_volume/PH_single/homer'
% num_levels = 1;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'mug_volume/PH_multi/mug'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'noisey_bunny_volume/PH_multi/noisey_bunny'
% num_levels = 4;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% % Obs.: Level 5 has geometric degenerate faces. So I'm computing
% % intersection only for the first 4
% 
% basename = 'octopus-300k_volume/PH_multi/octopus-300k'
% num_levels = 11;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% basename = 'alien_varap/PH_multi/alien'
% num_levels = 5;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'anchor_volume/PH_multi/anchor'
% num_levels = 6;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'arma_volumetric_arap/PH_multi/arma'
% num_levels = 6;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'bimba_varap/PH_multi/bimba'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'bunny-50k_varap/PH_multi/bunny'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'couplingdown_volume/PH_multi/couplingdown'
% num_levels = 3;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'disney1_varap/PH_multi/disney1'
% num_levels = 5;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'hand_varap/PH_multi/hand'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'horse_varap_25/PH_multi/horse'
% num_levels = 21;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% % Obs.: geometrically degenerate faces at level 22. So testing only the
% % first 21.
% 
% basename = 'lobster_varap/PH_multi/lobster'
% num_levels = 5;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'manhead_varap/PH_multi/manhead'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'maxplank_varap_50/PH_multi/maxplank'
% num_levels = 50;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'Model1_varap/PH_multi/Model1'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'Model3_varap/PH_multi/Model3'
% num_levels = 4;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'Model4_varap/PH_multi/Model4'
% num_levels = 4;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'Model7_volume/PH_multi/Model7'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');
% 
% basename = 'Model9_varap/PH_multi/Model9'
% num_levels = 7;
% [IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
% save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');

basename = 'pelvis_volume/PH_multi/pelvis'
num_levels = 7;
[IF_input,IF_levels] = progressive_hulls_intersections(basename,num_levels);
save(sprintf('%s_intersections.mat',basename),'IF_input','IF_levels');