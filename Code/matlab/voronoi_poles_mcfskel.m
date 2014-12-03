function [V_poles,F_poles] = voronoi_poles_mcfskel(V0,F0)
  % VORONOI_POLES_MCFSKEL
  % voronoi_poles_mcfskel(V0,F0,varagin)
  %
  % Runs on command line the filter 'Voronoi based MAT' from the starterm
  % application, which is the implentation of "Mean Curvature Skeleletons"
  % by Tagliasacchi et al. (SGP 2012).
  %
  % Inputs:
  %   V0  list of surface vertex positions of exterior mesh, # vertices by 3
  %   F0  list of surface face indices of exterior triangle mesh, # faces by 3
  % Outputs:
  %   V_poles  list of poles for each vertex of the initial mesh
  %   F_poles  (= F0)
  
  % get a temporary file name prefix
  prefix = tempname;
  
  % add extension.off to the name and write (V0,F0) to file
  off_filename = [prefix '.off'];
  writeOFF(off_filename,V0,F0);
  
  % set path to the skeletonization program
  path_to_mcfskel = '/Applications/Starlab.app/Contents/MacOS/starterm';
  
  % call starterm from command line
  command = [path_to_mcfskel ' ' '--filter="Voronoi based MAT"' ' ' off_filename];
  fprintf(command);
  [status, result] = system(command);
  
  % detele created OFF file
  delete(off_filename);
  
  % read pole file
  pole_file = fopen('poles.txt');
  V_poles = fscanf(pole_file,'%f',[3,inf])';
  fclose(pole_file);
  F_poles = F0;
  
  % delete pole file
  system('rm poles.txt');