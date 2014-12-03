function [J,Dg_Dp,Dp_Dv] = response_constraint_matrix_3D(F,Col,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,...
    EdgeCol,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2)
  % RESPONSE_CONSTRAINT_MATRIX_3D
  % response_constraint_matrix_3D(F,Col,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,...
  % EdgeCol,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2)
  %
  % Implements the matrix J as descibed in "Implicit Contact Handling for Deformable Objects" 
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
  % Output:
  %   J   (number of collisions)x(3*M) matrix with the coefficients of the
  %   constraint matrix J of Otaduy et al. 09 (section 3.2)
  %   Dg_Dp   (number of collisions)x(6*number of collisions) matrix 
  %   that has as entries collision normals
  %   Dp_Dv   (6*number of collisions)x(3*M) matrix relates velocities
  %   at collision points to velocities of mesh vertices
  
  % the matrix J will have num_collisions by 3*size(mesh),
  % to consider x, y and z velocities. We consider the unknowns to have
  % first all x velocities, then all y velocities and then all z velocities
  num_collisions = full(sum(sum(Col)));
  edge_num_collisions = full(sum(sum(EdgeCol)));
  
  J = sparse(num_collisions+edge_num_collisions,3*size(Col,1));
  
  Dg_Dp = sparse(num_collisions+edge_num_collisions,6*(num_collisions+edge_num_collisions));
  Dp_Dv = sparse(6*(num_collisions+edge_num_collisions),3*size(Col,1));
  
  [ColI,ColJ,~] = find(Col);
  [EdgeColI,EdgeColJ,~] = find(EdgeCol);
  
  % Beginning of defining matrix Dg_Dp that contains the normals
  I_Dg_Dp = zeros(1,6*(num_collisions+edge_num_collisions),1);
  J_Dg_Dp = zeros(1,6*(num_collisions+edge_num_collisions));
  Val_Dg_Dp = zeros(1,6*(num_collisions+edge_num_collisions));
  
  for k=1:num_collisions
      
      i = ColI(k); j = ColJ(k);
      
      I_Dg_Dp(6*(k-1)+1) = k;
      J_Dg_Dp(6*(k-1)+1) = 2*(k-1)+1;
      Val_Dg_Dp(6*(k-1)+1) = ColNX(i,j);
      
      I_Dg_Dp(6*(k-1)+2) = k;
      J_Dg_Dp(6*(k-1)+2) = 2*(k-1)+2;
      Val_Dg_Dp(6*(k-1)+2) = -ColNX(i,j);
      
      I_Dg_Dp(6*(k-1)+3) = k;
      J_Dg_Dp(6*(k-1)+3) = 2*(k-1)+1+2*num_collisions;
      Val_Dg_Dp(6*(k-1)+3) = ColNY(i,j);
      
      I_Dg_Dp(6*(k-1)+4) = k;
      J_Dg_Dp(6*(k-1)+4) = 2*(k-1)+2+2*num_collisions;
      Val_Dg_Dp(6*(k-1)+4) = -ColNY(i,j);
      
      I_Dg_Dp(6*(k-1)+5) = k;
      J_Dg_Dp(6*(k-1)+5) = 2*(k-1)+1+4*num_collisions;
      Val_Dg_Dp(6*(k-1)+5) = ColNZ(i,j);
      
      I_Dg_Dp(6*(k-1)+6) = k;
      J_Dg_Dp(6*(k-1)+6) = 2*(k-1)+2+4*num_collisions;
      Val_Dg_Dp(6*(k-1)+6) = -ColNZ(i,j);
      
  end
  
  for k=1:edge_num_collisions
      
      i = EdgeColI(k); j = EdgeColJ(k);
      
      I_Dg_Dp(6*num_collisions+6*(k-1)+1) = num_collisions+k;
      J_Dg_Dp(6*num_collisions+6*(k-1)+1) = num_collisions+2*(k-1)+1;
      Val_Dg_Dp(6*num_collisions+6*(k-1)+1) = EdgeColNX(i,j);
      
      I_Dg_Dp(6*num_collisions+6*(k-1)+2) = num_collisions+k;
      J_Dg_Dp(6*num_collisions+6*(k-1)+2) = num_collisions+2*(k-1)+2;
      Val_Dg_Dp(6*num_collisions+6*(k-1)+2) = -EdgeColNX(i,j);
      
      I_Dg_Dp(6*num_collisions+6*(k-1)+3) = num_collisions+k;
      J_Dg_Dp(6*num_collisions+6*(k-1)+3) = num_collisions+2*(k-1)+1+2*edge_num_collisions;
      Val_Dg_Dp(6*num_collisions+6*(k-1)+3) = EdgeColNY(i,j);
      
      I_Dg_Dp(6*num_collisions+6*(k-1)+4) = num_collisions+k;
      J_Dg_Dp(6*num_collisions+6*(k-1)+4) = num_collisions+2*(k-1)+2+2*edge_num_collisions;
      Val_Dg_Dp(6*num_collisions+6*(k-1)+4) = -EdgeColNY(i,j);
      
      I_Dg_Dp(6*num_collisions+6*(k-1)+5) = num_collisions+k;
      J_Dg_Dp(6*num_collisions+6*(k-1)+5) = num_collisions+2*(k-1)+1+4*edge_num_collisions;
      Val_Dg_Dp(6*num_collisions+6*(k-1)+5) = EdgeColNZ(i,j);
      
      I_Dg_Dp(6*num_collisions+6*(k-1)+6) = num_collisions+k;
      J_Dg_Dp(6*num_collisions+6*(k-1)+6) = num_collisions+2*(k-1)+2+4*edge_num_collisions;
      Val_Dg_Dp(6*num_collisions+6*(k-1)+6) = -EdgeColNZ(i,j);
      
  end

  Dg_Dp = sparse(I_Dg_Dp,J_Dg_Dp,Val_Dg_Dp,num_collisions+edge_num_collisions,6*(num_collisions+edge_num_collisions));
  
  % End of defining matrix Dg_Dp
  
  % Beginning of defining matrix Dp_Dv that contains barycentric
  % coordinates of the collision points w.r.t. edges
  % (this expresses velocities of colliding points in terms of velocities
  % of vertices of the mesh)
  
  I_Dp_Dv = zeros(1,12*(num_collisions+edge_num_collisions),1);
  J_Dp_Dv = zeros(1,12*(num_collisions+edge_num_collisions));
  Val_Dp_Dv = zeros(1,12*(num_collisions+edge_num_collisions));
  
  for k=1:num_collisions
      
      i = ColI(k); j = ColJ(k);
      j1 = F(j,1); j2 = F(j,2); j3 = F(j,3);
      u1 = ColU1(i,j); u2 = ColU2(i,j); u3 = ColU3(i,j);
      
      % first vertex of a colliding face
      I_Dp_Dv(12*(k-1)+1) = 2*(k-1)+1;
      J_Dp_Dv(12*(k-1)+1) = j1;
      Val_Dp_Dv(12*(k-1)+1) = u1;
      
      % second vertex of a colliding face
      I_Dp_Dv(12*(k-1)+2) = 2*(k-1)+1;
      J_Dp_Dv(12*(k-1)+2) = j2;
      Val_Dp_Dv(12*(k-1)+2) = u2;
      
      % third vertex of a colliding face
      I_Dp_Dv(12*(k-1)+3) = 2*(k-1)+1;
      J_Dp_Dv(12*(k-1)+3) = j3;
      Val_Dp_Dv(12*(k-1)+3) = u3;
      
      % for vertices the barycentric coordinate has to be one
      I_Dp_Dv(12*(k-1)+4) = 2*(k-1)+2;
      J_Dp_Dv(12*(k-1)+4) = i;
      Val_Dp_Dv(12*(k-1)+4) = 1.0;
      
      % Below it is all the same for y velocities
      I_Dp_Dv(12*(k-1)+5) = 2*(k-1)+1+2*num_collisions;
      J_Dp_Dv(12*(k-1)+5) = j1+size(Col,1);
      Val_Dp_Dv(12*(k-1)+5) = u1;
      
      I_Dp_Dv(12*(k-1)+6) = 2*(k-1)+1+2*num_collisions;
      J_Dp_Dv(12*(k-1)+6) = j2+size(Col,1);
      Val_Dp_Dv(12*(k-1)+6) = u2;
      
      I_Dp_Dv(12*(k-1)+7) = 2*(k-1)+1+2*num_collisions;
      J_Dp_Dv(12*(k-1)+7) = j3+size(Col,1);
      Val_Dp_Dv(12*(k-1)+7) = u3;
      
      I_Dp_Dv(12*(k-1)+8) = 2*(k-1)+2+2*num_collisions;
      J_Dp_Dv(12*(k-1)+8) = i+size(Col,1);
      Val_Dp_Dv(12*(k-1)+8) = 1.0;
      
      % Below it is all the same for z velocities
      I_Dp_Dv(12*(k-1)+9) = 2*(k-1)+1+4*num_collisions;
      J_Dp_Dv(12*(k-1)+9) = j1+2*size(Col,1);
      Val_Dp_Dv(12*(k-1)+9) = u1;
      
      I_Dp_Dv(12*(k-1)+10) = 2*(k-1)+1+4*num_collisions;
      J_Dp_Dv(12*(k-1)+10) = j2+2*size(Col,1);
      Val_Dp_Dv(12*(k-1)+10) = u2;
      
      I_Dp_Dv(12*(k-1)+11) = 2*(k-1)+1+4*num_collisions;
      J_Dp_Dv(12*(k-1)+11) = j3+2*size(Col,1);
      Val_Dp_Dv(12*(k-1)+11) = u3;
      
      I_Dp_Dv(12*(k-1)+12) = 2*(k-1)+2+4*num_collisions;
      J_Dp_Dv(12*(k-1)+12) = i+2*size(Col,1);
      Val_Dp_Dv(12*(k-1)+12) = 1.0;
      
  end
  
  % compute all edges
  EE_all = edges(F);
  
  for k=1:edge_num_collisions
      
      i = EdgeColI(k); j = EdgeColJ(k);
      i1 = EE_all(i,1); i2 = EE_all(i,2); 
      j1 = EE_all(j,1); j2 = EE_all(j,2); 
      u1 = EdgeColU1(i,j); u2 = EdgeColU2(i,j);
      
      % first vertex of the first colliding edge
      I_Dp_Dv(12*num_collisions+12*(k-1)+1) = 6*num_collisions+2*(k-1)+1;
      J_Dp_Dv(12*num_collisions+12*(k-1)+1) = i1;
      Val_Dp_Dv(12*num_collisions+12*(k-1)+1) = 1-u1;
      
      % second vertex of the first colliding edge
      I_Dp_Dv(12*num_collisions+12*(k-1)+2) = 6*num_collisions+2*(k-1)+1;
      J_Dp_Dv(12*num_collisions+12*(k-1)+2) = i2;
      Val_Dp_Dv(12*num_collisions+12*(k-1)+2) = u1;
      
      % first vertex of the second colliding edge
      I_Dp_Dv(12*num_collisions+12*(k-1)+3) = 6*num_collisions+2*(k-1)+2;
      J_Dp_Dv(12*num_collisions+12*(k-1)+3) = j1;
      Val_Dp_Dv(12*num_collisions+12*(k-1)+3) = 1-u2;
      
      % second vertex of the second colliding edge
      I_Dp_Dv(12*num_collisions+12*(k-1)+4) = 6*num_collisions+2*(k-1)+2;
      J_Dp_Dv(12*num_collisions+12*(k-1)+4) = j2;
      Val_Dp_Dv(12*num_collisions+12*(k-1)+4) = u2;
      
      % Below it is all the same for y velocities
      I_Dp_Dv(12*num_collisions+12*(k-1)+5) = 6*num_collisions+2*(k-1)+1+2*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+5) = i1+size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+5) = 1-u1;
      
      I_Dp_Dv(12*num_collisions+12*(k-1)+6) = 6*num_collisions+2*(k-1)+1+2*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+6) = i2+size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+6) = u1;
      
      I_Dp_Dv(12*num_collisions+12*(k-1)+7) = 6*num_collisions+2*(k-1)+2+2*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+7) = j1+size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+7) = 1-u2;
      
      I_Dp_Dv(12*num_collisions+12*(k-1)+8) = 6*num_collisions+2*(k-1)+2+2*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+8) = j2+size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+8) = u2;
      
      % Below it is all the same for z velocities
      I_Dp_Dv(12*num_collisions+12*(k-1)+9) = 6*num_collisions+2*(k-1)+1+4*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+9) = i1+2*size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+9) = 1-u1;
      
      I_Dp_Dv(12*num_collisions+12*(k-1)+10) = 6*num_collisions+2*(k-1)+1+4*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+10) = i2+2*size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+10) = u1;
      
      I_Dp_Dv(12*num_collisions+12*(k-1)+11) = 6*num_collisions+2*(k-1)+2+4*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+11) = j1+2*size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+11) = 1-u2;
      
      I_Dp_Dv(12*num_collisions+12*(k-1)+12) = 6*num_collisions+2*(k-1)+2+4*edge_num_collisions;
      J_Dp_Dv(12*num_collisions+12*(k-1)+12) = j2+2*size(Col,1);
      Val_Dp_Dv(12*num_collisions+12*(k-1)+12) = u2;
      
  end
  
  Dp_Dv = sparse(I_Dp_Dv,J_Dp_Dv,Val_Dp_Dv,6*(num_collisions+edge_num_collisions),3*size(Col,1));
  
  % End of defining matrix Dp_Dv
  
  % final matrix is the product of the matrices
  J = Dg_Dp*Dp_Dv;
  
  