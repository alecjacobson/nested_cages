function [SD,DV,DT,DF] = distance_on_conforming_tetmesh(V0,F0)

  % DISTANCE_ON_CONFORMING_TETMESH
  % [SD,DV,DT,DF] = distance_on_conforming_tetmesh(V0,F0)
  %
  % Calculate a tetmesh for a box that contains a surface (V0,F0) where
  % the vertices of (V0,F0) also belong to the tet-mesh. Additionally,
  % it computes signed distances w.r.t. (V0,F0) for the whole tetmesh.
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the input mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   input mesh
  % Output:
  %   (DV,DT,DF): tet-mesh that contains the input surface
  %   SD: signed distances of the vertices of this tetmesh w.r.t. input
  %   surface.
  
  % Compute tetmesh (tune parameters later if needed)
  
  [DV,DT,DF] = cdt(V0,F0,'UseBoundingBox',true,'BoundingBoxPush',1.5,'TetgenFlags','-q2.0 -a1e-4');
  
  % calculate (absoulte) distances of the vertices of the tet-mesh 
  % w.r.t. (V0,F0)
  D = sqrt(point_mesh_squared_distance(DV,V0,F0));
  
  % calculate winding number to determine signs of the distances
  W = winding_number(V0,F0,DV)/(4*pi); % this belongs to {0,1}

  % multiply signs by distances
  SD = (W*2-1).*D;
  % flipping seems the correct thing to do
  SD = -SD;