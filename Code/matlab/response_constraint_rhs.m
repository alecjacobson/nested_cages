function [g0] = response_constraint_rhs(E,Col,ColNX,ColNY,ColU,V0,delta_t)
  % RESPONSE_CONSTRAINT_RHS
  % response_constraint_rhs(E,Col,ColNX,ColNY,ColU,V0,delta_t)
  %
  % Implements the vector g0 as descibed in "Implicit Contact Handling for Deformable Objects" 
  % by Otaduy et al. 2009 (section 3.2)
  %
  % Input:
  %   E  2 x N list of vertex indices that form the edges of the mesh 
  %   Col   MxN matrix such that Col(i,j) == 1 iff there is a collision between
  %     vertex i and edges j 
  %   ColN1   MxN matrix such that ColN1(i,j) is the x-coordinate
  %     of the collision normal of vertex i collides with edge j. 
  %   ColN2   MxN matrix such that ColN2(i,j) is the y-coordinate
  %     of the collision normal of vertex i collides with edge j. 
  %   ColU   MxN matrix such that ColU(i,j) is the fraction of edge
  %      j where vertex i collides. 
  %   V0   Nx2 list of vertex positions at time t0
  %   delta_t   time step
  % Output:
  %   g0   (number of collisions)x1 vector entries of the
  %   vector g0 of Otaduy et al. 09 (section 3.2)
   
  
  % indicies of the collisions
  [ColI,ColJ,~] = find(Col);
  % collision normals (at time t1)
  Normals = [ColNX(ColI+size(V0,1)*(ColJ-1)) ColNY(ColI+size(V0,1)*(ColJ-1))];
  
%   % initial normals
%   Edges_0 = V0(E(:,1),:)-V0(E(:,2),:);
%   Normals_0 = [-Edges_0(:,2) Edges_0(:,1)];
%   Normals = Normals_0(E(ColJ,1),:);
%   Normals = normalizerow(Normals);
  
%   % final normals
%   Edges_1 = V1(E(:,1),:)-V1(E(:,2),:);
%   Normals_1 = [-Edges_1(:,2) Edges_1(:,1)];
%   Normals = Normals_1(E(ColJ,1),:);
%   Normals = normalizerow(Normals);
  
  % points on colliding edges at time t0
  Pa_0 = [[ColU(ColI+size(V0,1)*(ColJ-1)) ColU(ColI+size(V0,1)*(ColJ-1))].*V0(E(ColJ,2),:)...
      + [1-ColU(ColI+size(V0,1)*(ColJ-1)) 1-ColU(ColI+size(V0,1)*(ColJ-1))].*V0(E(ColJ,1),:)];
  % colliding vertices at time t0
  Pb_0 = V0(ColI,:);
  
%   % CCW test to know if it necessary to flip normals
%   ccw_test = zeros(size(ColI,1),1);
%   for k=1:size(ColI,1)
%       ccw_test(k) = det([V0(ColI(k),1) V0(E(ColJ(k),1),1) V0(E(ColJ(k),2),1);...
%                            V0(ColI(k),2) V0(E(ColJ(k),1),2) V0(E(ColJ(k),2),2);1 1 1]);
%   end
  
  % final constraint rhs (+plus flipping normals when necessary) - old way
%   g0 = (-1/delta_t)*sum(Normals.*(Pa_0-Pb_0),2).*ccw_test;
  % final constraint rhs - right way
    g0 = (-1/delta_t)*sum(Normals.*(Pa_0-Pb_0),2)+0.001/delta_t;
%   g0 = (-1/delta_t)*sum(Normals.*(Pa_0-Pb_0),2)+0.0001/delta_t;