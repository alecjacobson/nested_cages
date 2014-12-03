function [origin,distances,dx] = readSDF(sdf_file_name)

  % READSDF signed distances from a .SDF file
  %
  % [origin,distances,dx] = readSDF(sdf_file_name)
  %
  % Input:
  %  sdf_file_name  path to .sdf file
  %
  % Output:
  %   origin  3D position of the of the leftmost, bottom most corner of the
  %   bounding box
  %   distances  num_cells_x x num_cells_y x num_cells_z by 1 vector with
  %   signed distances for each cell
  %   dx  3D voxel dimension size
    
    sdf_file_handle = fopen(sdf_file_name);
    dimensions = fscanf(sdf_file_handle,'%d %d %d',3);
    nx = dimensions(1);
    ny = dimensions(2);
    nz = dimensions(3);
    origin = fscanf(sdf_file_handle,'%f %f %f',3);
    dx = fscanf(sdf_file_handle,'%f',1);
    
    distances = zeros(nx,ny,nz);
    
    for k =1:nz
        for j=1:ny
            for i=1:nx
                distances(i,j,k) = fscanf(sdf_file_handle,'%f',1);
            end
        end
    end