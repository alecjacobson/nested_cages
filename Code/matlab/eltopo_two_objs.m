function V_final = eltopo_two_objs(V0,F0,V1,N,varargin)
% ELTOPO_TWO_OBJS
% eltopo_two_objs(V0,F0,V1,N,varargin)
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

% get a temporary file name prefix0 and define OBJ file
  prefix0 = tempname;
  obj_filename0 = [prefix0 '.obj'];
  obj_filename0 = '/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/meshes/V0_prob1.obj';
  writeOBJ(obj_filename0,V0,F0);
% get a temporary file name prefix1
  prefix1 = tempname;
  obj_filename1 = [prefix1 '.obj'];
  obj_filename1 = '/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/meshes/V1_prob1.obj';
  writeOBJ(obj_filename1,V1,F0);
 % get a temporary file name prefix_out
  prefix_out = tempname;
  obj_filename_out = [prefix_out '.obj'];
  
  command = ['/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/collide_eltopo_two_objs'...
      ' ' obj_filename0 ' ' obj_filename1 ' ' sprintf('%d',N) ' ' obj_filename_out]
  [status, result] = system(command);
  
  disp('ended simulation. Now read obj')
  
  if status~=0
    error(result)
  end
  
  [V_final,~] = readOBJ(obj_filename_out);
  
  % delete all files
%   delete(obj_filename0);
%   delete(obj_filename1);
  delete(obj_filename_out);



