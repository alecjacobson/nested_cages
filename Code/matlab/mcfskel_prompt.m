function [V,F] = mcfskel_prompt(V0,F0)
  % MCFSKEL_PROMPT
  % [V,F] = mcfskel_prompt(V0,F0)
  %
  % Calls on the command line the program available as 
  % implementation of the paper "Mean Curvature Skeleletons"
  % by Tagliasacchi et al. (SGP 2012).
  % The program performs only one contraction of the mesh towards
  % the medial axis.
  %
  % Inputs:
  %   V0  list of surface vertex positions of exterior mesh, # vertices by 3
  %   F0  list of surface face indices of exterior triangle mesh, # faces by 3
  % Outputs:
  %   V  list of surface vertices
  %   F  list of faces (different from the input one)
    
  % get a temporary file name prefix
  prefix = tempname;
  
  % add extension.off to the name and write (V0,F0) to file
  off_filename = [prefix '.off'];
  off_filename_ckel = [prefix '_ckel.off'];
  writeOFF(off_filename,V0,F0);
  
  % set path to the skeletonization program
  path_to_mcfskel = '/Applications/Starlab.app/Contents/MacOS/starterm';
  
  % I'm not being able to call mcfskel from Matlab. 
  % So I just print the command, copy it and paste on the Terminal.
  % This is the only way I found
  
  % call MCFskel and overwrite the result
  command = [path_to_mcfskel ' ' '--filter="MCF Skeletonization"' ' ' off_filename ' ' '--save-overwrite']
%   [status, result] = system(command);
%   result

  input('');
    
  % read result
  [V,F] = readOFF(off_filename_ckel);
  
  % detele created OFF file
  delete(off_filename);
  delete(off_filename_ckel);