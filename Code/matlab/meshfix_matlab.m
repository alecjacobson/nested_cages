function [V,F] = meshfix_matlab(V0,F0)
  % MESHFIX_MATLAB
  % meshfix_matlab(V0,F0)
  %
  % Given a triangle mesh (V0,F0) calls MeshFix binary to fix it.
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the mesh
  %     
  % Output:
  %   (V,F)   fixed mesh (output of MeshFix)
  
  prefix = tempname;
  off_filename = [prefix '.off'];
  writeOFF(off_filename,V0,F0);
  
  command = ['/usr/bin/meshfix' ' ' off_filename];
  fprintf(command);
  [status, result] = system(command);
  
  output_filename = [prefix '_fixed.off']
  [V,F] = readOFF(output_filename);
  
  delete(off_filename);