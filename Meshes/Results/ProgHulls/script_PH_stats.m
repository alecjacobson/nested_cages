basename = 'alien_varap/PH_multi/alien';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'anchor_volume/PH_multi/anchor';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'arma_volumetric_arap/PH_multi/arma';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'bimba_varap/PH_multi/bimba';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'bunny-50k_varap/PH_multi/bunny';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'couplingdown_volume/PH_multi/couplingdown';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'disney1_varap/PH_multi/disney1';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'fert_volume/PH_multi/fert';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'gallop_arap/PH_single/horse_001';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'gargo_volume/PH_multi/gargo';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'hand_varap/PH_multi/hand';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'handles_volume/PH_multi/handles';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'homer_volume/PH_multi/homer';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'horse_varap_25/PH_multi/horse';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'lobster_varap/PH_multi/lobster';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'manhead_varap/PH_multi/manhead';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'maxplank_varap_50/PH_multi/maxplank';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'Model1_varap/PH_multi/Model1';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'Model3_varap/PH_multi/Model3';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'Model4_varap/PH_multi/Model4';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'Model7_volume/PH_multi/Model7';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'Model9_varap/PH_multi/Model9';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'mug_volume/PH_multi/mug';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'noisey_bunny_volume/PH_multi/noisey_bunny';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'octopus-300k_volume/PH_multi/octopus-300k';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

basename = 'pelvis_volume/PH_multi/pelvis';
load(sprintf('%s_intersections.mat',basename));
num_levels = size(diag(IF_levels),1);
self_int = sum(diag(IF_levels)>0);
input_int = sum(IF_input>0);
inter_levels_int = sum(sum(triu(IF_levels,1)>0));
inter_levels_total = (num_levels-1)*(num_levels)/2;
message = sprintf('basename = %s \ninput: %d/%d - %.2f%% \nself: %d/%d - %.2f%% \nlevels: %d/%d - %.2f%%',basename,...
    input_int,num_levels, 100*(input_int/num_levels),...
    self_int,num_levels, 100*(self_int/num_levels),...
    inter_levels_int,inter_levels_total,100*(inter_levels_int/inter_levels_total));
disp(message)

% .                       anchor_volume           horse_varap_25          
% ..                      animal_arap             lobster_varap           
% .DS_Store               arma_volumetric_arap    manhead_varap           
% Model1_varap            bimba_varap             maxplank_varap_50       
% Model3_varap            bunny-50k_varap         mug_volume              
% Model4_005              couplingdown_volume     noisey_bunny_volume     
% Model4_05               disney1_varap           octopus-300k_volume     
% Model4_varap            fert_volume             pelvis_volume           
% Model7_volume           gallop_arap             script_PH.m             
% Model9_005              gargo_volume            script_PH_stats.m       
% Model9_varap            hand_varap              script_PH_stats.m~      
% Model9_volume           handles_volume          script_intersections.m  
% alien_varap             homer_volume     