function V_final = eltopo_two_bins(V0,F0,V1,N,varargin)
% ERROR: THIS THING DOESN'T WORK!!! DON'T USE IT!
% ELTOPO_TWO_BINS
% eltopo_two_bins(V0,F0,V1,N,varargin)
%
% Given two meshes with the same connectivity (V0,F0) and (V1,F0) and a
% a number N that determines that the first N vertices have infinite mass,
% do collision detection and response using eltopo library from command
% line.
%
% Input:
%   V0 (#vertices)x3 list of mesh vertex positions of the mesh at the
%   beginning of the time step
%   F0 (#faces)x3 list of vertex indices that form each face of both
%   meshes
%   V1 (#vertices)x3 list of mesh vertex positions of the mesh at the
%   beginning of the time step
%   N: number of vertices that are constrained from one step to the other
%   (infinite mass points)
%   Optional:
%     
% Output:
%   V_final (#vertices)x3 list of mesh vertex positions of the
%   mesh after collision detection and response

% get a temporary file name prefix0 and define BIN file
  prefix0 = tempname;
  bin_filename0 = [prefix0 '.bin'];
  bin_filename0 = '/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/meshes/V0_prob1.bin';
  fid_0 = fopen(bin_filename0, 'a+');
%   size(V0,1)
%   fwrite(fid_0, size(V0,1), 'integer*4');
  F0
%   fwrite(fid_0, single(V0), 'float');
%   fwrite(fid_0, size(F0,1), 'integer*4');
  fwrite(fid_0, F0+100,'float');
  fclose(fid_0);
% get a temporary file name prefix1
  prefix1 = tempname;
  bin_filename1 = [prefix1 '.bin'];
  bin_filename1 = '/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/meshes/V1_prob1.bin';
  fid_1 = fopen(bin_filename1, 'w+');
  fwrite(fid_1, size(V1,1), 'integer*4');
  fwrite(fid_1, single(V1), 'float');
  fwrite(fid_1, size(F0,1), 'integer*4');
  fwrite(fid_1, F0, 'integer*4');
  fclose(fid_1);
  input('');
 % get a temporary file name prefix_out
  prefix_out = tempname;
  obj_filename_out = [prefix_out '.obj'];
  
  command = ['/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/collide_eltopo_two_objs'...
      ' ' bin_filename0 ' ' bin_filename1 ' ' sprintf('%d',N) ' ' obj_filename_out]
  [status, result] = system(command);
  if status~=0
    error(result)
  end
  
  [V_final,~] = readOBJ(obj_filename_out);
  
  % delete all files
  delete(bin_filename0);
  delete(bin_filename1);
  delete(obj_filename_out);