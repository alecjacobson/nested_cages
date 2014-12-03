function [g0] = response_constraint_rhs_3D(F,Col,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,...
    EdgeCol,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2,EdgeColT,V0,delta_t)
  % RESPONSE_CONSTRAINT_RHS_3D
  % response_constraint_rhs_3D(F,Col,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,...
  %  EdgeCol,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2,V0,delta_t)
  %
  % Implements the vector g0 as descibed in "Implicit Contact Handling for Deformable Objects" 
  % by Otaduy et al. 2009 (section 3.2)
  %
  % Input:
  %   F  N x 3 list of vertex indices that form the edges of the mesh 
  %   Col   MxN matrix such that Col(i,j) == 1 iff there is a collision between
  %     vertex i and face j 
  %   ColNX   MxN matrix such that ColNX(i,j) is the x-coordinate
  %     of the collision normal of vertex i collides with face j. 
  %   ColNY   MxN matrix such that ColNY(i,j) is the y-coordinate
  %     of the collision normal of vertex i collides with face j. 
  %   ColNZ   MxN matrix such that ColNZ(i,j) is the z-coordinate
  %     of the collision normal of vertex i collides with face j. 
  %   ColU1   MxN matrix such that ColU1(i,j) is the barycentric
  %   coordinate of the collision w.r.t. the 1st vertex
  %   ColU2   MxN matrix such that ColU2(i,j) is the barycentric
  %   coordinate of the collision w.r.t. the 1st vertex 
  %   ColU3   MxN matrix such that ColU3(i,j) is the barycentric
  %   coordinate of the collision w.r.t. the 1st vertex 
  %   EdgeCol   (#edges)x(#edges) matrix such that Col(i,j) == 1 
  %     iff there is a collision between edge i and edge j
  %   EdgeColNX   (#edges)x(#edges) matrix such that ColNX(i,j) 
  %   is the x-coordinate of the collision normal of edge 
  %   i collides with edge j. 
  %   EdgeColNY   (#edges)x(#edges) matrix such that ColNY(i,j) 
  %   is the y-coordinate of the collision normal of edge 
  %   i collides with edge j.
  %   EdgeColNZ   (#edges)x(#edges) matrix such that ColNZ(i,j) 
  %   is the z-coordinate of the collision normal of edge 
  %   i collides with edge j.
  %   EdgeColU1   (#edges)x(#edges) matrix such that ColU1(i,j) 
  %   is the fraction of edge i where it collides with edge j.
  %   EdgeColU2   (#edges)x(#edges) matrix such that ColU1(i,j) 
  %   is the fraction of edge j where it collides with edge i.
  %   EdgeColT   (#edges)x(#edges) matrix such that EdgeColT(i,j) 
  %   is the time instant when edge i collides with edge j.
  %   V0   Nx3 list of vertex positions at time t0
  %   V1   Nx3 list of vertex positions at time t0+delta_t
  %   delta_t   time step
  % Output:
  %   g0   (number of collisions)x1 vector entries of the
  %   vector g0 of Otaduy et al. 09 (section 3.2)
   
  
  % indicies of the collisions
  [ColI,ColJ,~] = find(Col);
  
%   % normals at time zero corresponding to collision normals 
%   Normals = normals(V0,F);
% %   Normals = normalizerow(Normals);
%   Normals = Normals(ColJ,:);
    Normals = zeros(size(ColI,1),3);

  for k=1:size(ColI,1)
      
      i = ColI(k);
      j = ColJ(k);
      
      Normals(k,:) = [ColNX(i,j) ColNY(i,j) ColNZ(i,j)];
      
  end
  
  % points on colliding faces at time t0
  Pa_0 = [[ColU1(ColI+size(V0,1)*(ColJ-1)) ColU1(ColI+size(V0,1)*(ColJ-1)) ColU1(ColI+size(V0,1)*(ColJ-1))].*V0(F(ColJ,1),:)...
      + [ColU2(ColI+size(V0,1)*(ColJ-1)) ColU2(ColI+size(V0,1)*(ColJ-1)) ColU2(ColI+size(V0,1)*(ColJ-1))].*V0(F(ColJ,2),:) ...
      + [ColU3(ColI+size(V0,1)*(ColJ-1)) ColU3(ColI+size(V0,1)*(ColJ-1)) ColU3(ColI+size(V0,1)*(ColJ-1))].*V0(F(ColJ,3),:)];
  % colliding vertices at time t0
  Pb_0 = V0(ColI,:);
  
  % final constraint rhs
  g0 = (-1/delta_t)*sum(Normals.*(Pa_0-Pb_0),2);

  % Now do all the same for edge-edge collisions
  % indicies of the collisions
  [EdgeColI,EdgeColJ,~] = find(EdgeCol);
  Normals = zeros(size(EdgeColI,1),3);
  
  % compute all edges
  EE_all = edges(F);
  
  for k=1:size(EdgeColI,1)
      i = EdgeColI(k);
      j = EdgeColJ(k);
      
      Normals(k,:) = [EdgeColNX(i,j) EdgeColNY(i,j) EdgeColNZ(i,j)];

  end

  I1 = EE_all(EdgeColI,1);
  I2 = EE_all(EdgeColI,2);
  J1 = EE_all(EdgeColJ,1);
  J2 = EE_all(EdgeColJ,2);
  % first point on colliding edges at time t0
  Pa_0 = [1-EdgeColU1(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) 1-EdgeColU1(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) 1-EdgeColU1(EdgeColI+size(EE_all,1)*(EdgeColJ-1))].*V0(I1,:)...
      + [EdgeColU1(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) EdgeColU1(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) EdgeColU1(EdgeColI+size(EE_all,1)*(EdgeColJ-1))].*V0(I2,:);
  % second point on colliding edges at time t0
  Pb_0 = [1-EdgeColU2(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) 1-EdgeColU2(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) 1-EdgeColU2(EdgeColI+size(EE_all,1)*(EdgeColJ-1))].*V0(J1,:)...
      + [EdgeColU2(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) EdgeColU2(EdgeColI+size(EE_all,1)*(EdgeColJ-1)) EdgeColU2(EdgeColI+size(EE_all,1)*(EdgeColJ-1))].*V0(J2,:);
  
  % final constraint rhs (+plus flipping normals when necessary)
  g0 = [g0; (-1/delta_t)*sum(Normals.*(Pa_0-Pb_0),2)]+1e-6/delta_t;  