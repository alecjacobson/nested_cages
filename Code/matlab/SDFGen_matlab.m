function [origin,distances,dx] = SDFGen_matlab(V,F,dx,padding)

  % SDFGen_matlab calls SDFGen to calculate the signed distance filed
  % with respect to given mesh.
  %
  % [origin,distances,dx] = SDFGen_matlab(V,F,dx,padding)
  %
  % Inputs:
  %   V  list of surface vertex positions of exterior mesh, # vertices by 3
  %   F  list of surface face indices of exterior triangle mesh, # faces by 3
  %   dx  3D voxel dimension size
  %   padding  number of extra cells out of the bounding box of the mesh
  % Outputs:
  %   origin  3D position of the of the leftmost, bottom most corner of the
  %   bounding box
  %   distances  num_cells_x x num_cells_y x num_cells_z by 1 vector with
  %   signed distances for each cell
  %   dx  3D voxel dimension size
  
  
  % Check if input variables are in cache (from cache_test.m)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Check for cached result, do NOT edit variables until cache is checked,
  % your function code comes later. See below
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get a list of current variables in this scope, this is the input "state"
  variables = who;
  % get a temporary file's name
  tmpf = [tempname('.') '.mat'];
  % save the "state" to file, so we can get a md5 checksum
  save(tmpf,'-regexp',sprintf('^%s$|',variables{:}),'-ascii');
  % get md5 checksum on input "state", we append .cache.mat to the check sum
  % because we'll use the checksum as the cache file name
  [s,cache_name] = system(['md5 -r ' tmpf ' | awk ''{printf "."$1".cache.mat"}''']);
  % clean up
  delete(tmpf);
  clear s tmpf variables;

  % If the checksum cache file exists then we've seen this input "state"
  % before, load in cached output "state"
  if(exist(cache_name,'file'))
    fprintf('Cache exists. Using cache...\n');
    % use cache
    load(cache_name);
  % Otherwise this is the first time we've seen this input "state", so we
  % execute the function as usual and save the output "state" to the cache 
  else
    fprintf('First time. Creating cache...\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Your function code goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get a temporary file name prefix
    prefix = tempname;
    
    % SDFGen accepts only OBJ files for triangle meshes
    obj_filename = [prefix '.obj'];
    writeOBJ(obj_filename,V,F);
    sdf_filename = [prefix '.sdf'];
    
    flags = sprintf('%f %d ',dx,padding);
    
    path_to_sdfgen = 'usr/bin/SDFGen';
    % call SDFGen
    command = [path_to_sdfgen ' ' obj_filename ' ' flags];
    
    [status, result] = system(command);
    if status~=0
        error(result)
    end
    
    [origin,distances,dx] = readSDF(sdf_filename);
    
    delete(obj_filename);
    delete(sdf_filename);
    
        % get list of variables present in this scope at finish of function code,
    % this is the output "state"
    variables = who;
    % save output "state" to file, using md5 checksum cache file name
    save(cache_name,'-regexp',sprintf('^%s$|',variables{:}));
  end