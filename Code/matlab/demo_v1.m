% demo_v1: loads a mesh, runs coarsening + shrinking, concatenates meshes
% and plot both of them togetther.

% read mesh
[V0 F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/arma.off');
%[V0 F0] = readOFF('/Users/Leo/PHD_Work/Cage_Generation_2013/models/rampant-R13k.off');

% run coarsening + shrinking 
[V_coarse F_coarse V_shrink F_shrink] = cage_generation_v1(V0,F0,300);

% concatenate results in a single mesh (only for visualization purposes)
V_all = [V_coarse;V_shrink];
F_all = [F_coarse;F_shrink+size(V_coarse,1)];

% visualize meshes (turn 'fill' option off and 'wireframe' on to see both)
meshplot(V_all,F_all);
