function [J,Dg_Dp,Dp_Dv] = response_constraint_matrix(E,Col,ColNX,ColNY,ColU)
  % RESPONSE_CONSTRAINT_MATRIX
  % response_constraint_matrix(E,Col,ColNX,ColNY,ColU)
  %
  % Implements the matrix J as descibed in "Implicit Contact Handling for Deformable Objects" 
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
  % Output:
  %   J   (number of collisions)x(2*M) matrix with the coefficients of the
  %   constraint matrix J of Otaduy et al. 09 (section 3.2)
  %   Dg_Dp   (number of collisions)x(4*number of collisions) matrix 
  %   that has as entries collision normals
  %   Dp_Dv   (4*number of collisions)x(2*M) matrix relates velocities
  %   at collision points to velocities of mesh vertices
  
  % the matrix J will have num_collisions by 2*size(mesh),
  % to consider x and y velocities. We consider the unknowns to have
  % first all x velocities and then all y velocities
  num_collisions = full(sum(sum(Col)));
  J = sparse(num_collisions,2*size(Col,1));
  
  Dg_Dp = sparse(num_collisions,4*num_collisions);
  Dp_Dv = sparse(4*num_collisions,2*size(Col,1));
  
  [ColI,ColJ,~] = find(Col);
  
  
  % Beginning of defining matrix Dg_Dp that contains the normals
  I_Dg_Dp = zeros(1,4*num_collisions,1);
  J_Dg_Dp = zeros(1,4*num_collisions);
  Val_Dg_Dp = zeros(1,4*num_collisions);
  
  for k=1:num_collisions
      
      i = ColI(k); j = ColJ(k);
      
      I_Dg_Dp(4*(k-1)+1) = k;
      J_Dg_Dp(4*(k-1)+1) = 2*(k-1)+1;
      Val_Dg_Dp(4*(k-1)+1) = ColNX(i,j);
      
      I_Dg_Dp(4*(k-1)+2) = k;
      J_Dg_Dp(4*(k-1)+2) = 2*(k-1)+2;
      Val_Dg_Dp(4*(k-1)+2) = -ColNX(i,j);
      
      I_Dg_Dp(4*(k-1)+3) = k;
      J_Dg_Dp(4*(k-1)+3) = 2*(k-1)+1+2*num_collisions;
      Val_Dg_Dp(4*(k-1)+3) = ColNY(i,j);
      
      I_Dg_Dp(4*(k-1)+4) = k;
      J_Dg_Dp(4*(k-1)+4) = 2*(k-1)+2+2*num_collisions;
      Val_Dg_Dp(4*(k-1)+4) = -ColNY(i,j);
      
  end

  Dg_Dp = sparse(I_Dg_Dp,J_Dg_Dp,Val_Dg_Dp,num_collisions,4*num_collisions);
  
  % End of defining matrix Dg_Dp
  
  % Beginning of defining matrix Dp_Dv that contains barycentric
  % coordinates of the collision points w.r.t. edges
  % (this expresses velocities of colliding points in terms of velocities
  % of vertices of the mesh)
  
  I_Dp_Dv = zeros(1,6*num_collisions,1);
  J_Dp_Dv = zeros(1,6*num_collisions);
  Val_Dp_Dv = zeros(1,6*num_collisions);
  
  for k=1:num_collisions
      
      i = ColI(k); j = ColJ(k);
      j1 = E(j,1); j2 = E(j,2);
      u = ColU(i,j);
      
      % first vertex of a colliding edge
      I_Dp_Dv(6*(k-1)+1) = 2*(k-1)+1;
      J_Dp_Dv(6*(k-1)+1) = j1;
      %Val_Dp_Dv(6*(k-1)+1) = u;
      Val_Dp_Dv(6*(k-1)+1) = (1-u);
      
      % second vertex of a colliding edge
      I_Dp_Dv(6*(k-1)+2) = 2*(k-1)+1;
      J_Dp_Dv(6*(k-1)+2) = j2;
      %Val_Dp_Dv(6*(k-1)+2) = (1-u);
      Val_Dp_Dv(6*(k-1)+2) = u;
      
      % for vertices the barycentric coordinate has to be one
      I_Dp_Dv(6*(k-1)+3) = 2*(k-1)+2;
      J_Dp_Dv(6*(k-1)+3) = i;
      Val_Dp_Dv(6*(k-1)+3) = 1.0;
      
      % Below it is all the same for y velocities
      I_Dp_Dv(6*(k-1)+4) = 2*(k-1)+1+2*num_collisions;
      J_Dp_Dv(6*(k-1)+4) = j1+size(Col,1);
      %Val_Dp_Dv(6*(k-1)+4) = u;
      Val_Dp_Dv(6*(k-1)+4) = (1-u);
      
      I_Dp_Dv(6*(k-1)+5) = 2*(k-1)+1+2*num_collisions;
      J_Dp_Dv(6*(k-1)+5) = j2+size(Col,1);
      %Val_Dp_Dv(6*(k-1)+5) = (1-u);
      Val_Dp_Dv(6*(k-1)+5) = u;
      
      I_Dp_Dv(6*(k-1)+6) = 2*(k-1)+2+2*num_collisions;
      J_Dp_Dv(6*(k-1)+6) = i+size(Col,1);
      Val_Dp_Dv(6*(k-1)+6) = 1.0;
      
  end
  
  Dp_Dv = sparse(I_Dp_Dv,J_Dp_Dv,Val_Dp_Dv,4*num_collisions,2*size(Col,1));
  
  % End of defining matrix Dp_Dv
  
%   full(Dg_Dp)
%   full(Dp_Dv)
  
  % final matrix is the product of the matrices
  J = Dg_Dp*Dp_Dv;  
  