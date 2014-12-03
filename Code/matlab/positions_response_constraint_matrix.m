function [N,B,A] = positions_response_constraint_matrix(E,Col,ColNX,ColNY,ColU,V_fine,E_fine)
  % POSITIONS_RESPONSE_CONSTRAINT_MATRIX
  % positions_response_constraint_matrix(E,Col,ColNX,ColNY,ColU)
  %
  % Implements constraint matrices N^{un} and N^{c} as decribed in my
  % report.
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
  %   (V_fine, E_fine) fine mesh
  % Output:
  %   N   (number of collisions)x(3*(number of collisions)) 
  %   matrix with entries being the entries of the collision normals
  %   B   (3*(number of collisions))x(3*num.vert. coarse mesh) barycentric
  %   coordinates of the vertices on the coarse mesh
  %   A   (3*(number of collisions))x(3*num.vert. fine mesh) barycentric
  %   coordinates of the vertices on the fine mesh
  
  num_collisions = full(sum(sum(Col)));
  
  % number of vertices in the fine (constrained) mesh
  nc = size(V_fine,1);
  % number of vertices in the coarse (unconstrained) mesh
  nun = size(Col,1)-nc;
    
  [ColI,ColJ,~] = find(Col);
  
  I_N = ones(1,2*num_collisions);
  J_N = ones(1,2*num_collisions);
  Val_N = zeros(1,2*num_collisions);
  
  I_B = ones(1,6*num_collisions);
  J_B = ones(1,6*num_collisions);
  Val_B = zeros(1,6*num_collisions);
  
  I_A = ones(1,4*num_collisions);
  J_A = ones(1,4*num_collisions);
  Val_A = zeros(1,4*num_collisions);
  
    for k=1:num_collisions
      
      i = ColI(k); j = ColJ(k);
      
      % case 1: vertex i belongs to fine mesh, face j to the coarse mesh
      if i<=nc
          
          j1 = E(j,1)-nc; j2 = E(j,2)-nc;

          % N matrix entries

          I_N(2*(k-1)+1) = k;
          J_N(2*(k-1)+1) = 2*(k-1)+1;
          Val_N(2*(k-1)+1) = ColNX(i,j);

          I_N(2*(k-1)+2) = k;
          J_N(2*(k-1)+2) = 2*(k-1)+2;
          Val_N(2*(k-1)+2) = ColNY(i,j);

          % B matrix entries

          I_B(6*(k-1)+1) = 2*(k-1)+1;
          J_B(6*(k-1)+1) = j1;
          Val_B(6*(k-1)+1) = 1-ColU(i,j);

          I_B(6*(k-1)+2) = 2*(k-1)+1;
          J_B(6*(k-1)+2) = j2;
          Val_B(6*(k-1)+2) = ColU(i,j);

          I_B(6*(k-1)+3) = 2*(k-1)+2;
          J_B(6*(k-1)+3) = j1+nun;
          Val_B(6*(k-1)+3) = 1-ColU(i,j);

          I_B(6*(k-1)+4) = 2*(k-1)+2;
          J_B(6*(k-1)+4) = j2+nun;
          Val_B(6*(k-1)+4) = ColU(i,j);

          % A matrix entries

          I_A(4*(k-1)+1) = 2*(k-1)+1;
          J_A(4*(k-1)+1) = i;
          Val_A(4*(k-1)+1) = 1.0;

          I_A(4*(k-1)+2) = 2*(k-1)+2;
          J_A(4*(k-1)+2) = i+nc;
          Val_A(4*(k-1)+2) = 1.0;

          
      % case 2: both vertex i and face j belong to coarse mesh    
      elseif (i>nc && j>size(E_fine,1))
          
          i = i-nc;
          j1 = E(j,1)-nc; j2 = E(j,2)-nc;
          
          % N matrix entries

          I_N(2*(k-1)+1) = k;
          J_N(2*(k-1)+1) = 2*(k-1)+1;
          Val_N(2*(k-1)+1) = ColNX(i+nc,j);

          I_N(2*(k-1)+2) = k;
          J_N(2*(k-1)+2) = 2*(k-1)+2;
          Val_N(2*(k-1)+2) = ColNY(i+nc,j);
          
          % B matrix entries

          I_B(6*(k-1)+1) = 2*(k-1)+1;
          J_B(6*(k-1)+1) = i;
          Val_B(6*(k-1)+1) = 1.0;

          I_B(6*(k-1)+2) = 2*(k-1)+2;
          J_B(6*(k-1)+2) = i+nun;
          Val_B(6*(k-1)+2) = 1.0;
          
          I_B(6*(k-1)+3) = 2*(k-1)+1;
          J_B(6*(k-1)+3) = j1;
          Val_B(6*(k-1)+3) = -(1-ColU(i+nc,j));

          I_B(6*(k-1)+4) = 2*(k-1)+1;
          J_B(6*(k-1)+4) = j2;
          Val_B(6*(k-1)+4) = -(ColU(i+nc,j));

          I_B(6*(k-1)+5) = 2*(k-1)+2;
          J_B(6*(k-1)+5) = j1+nun;
          Val_B(6*(k-1)+5) = -(1-ColU(i+nc,j));

          I_B(6*(k-1)+6) = 2*(k-1)+2;
          J_B(6*(k-1)+6) = j2+nun;
          Val_B(6*(k-1)+6) = -ColU(i+nc,j);

      
      % case 3: face j belongs to fine mesh, vertex i belongs to coarse
      else
          
          i = i-nc;
          j1 = E(j,1); j2 = E(j,2);
          
          % N matrix entries

          I_N(2*(k-1)+1) = k;
          J_N(2*(k-1)+1) = 2*(k-1)+1;
          Val_N(2*(k-1)+1) = ColNX(i+nc,j);

          I_N(2*(k-1)+2) = k;
          J_N(2*(k-1)+2) = 2*(k-1)+2;
          Val_N(2*(k-1)+2) = ColNY(i+nc,j);
          
          % B matrix entries

          I_B(6*(k-1)+1) = 2*(k-1)+1;
          J_B(6*(k-1)+1) = i;
          Val_B(6*(k-1)+1) = 1.0;

          I_B(6*(k-1)+2) = 2*(k-1)+2;
          J_B(6*(k-1)+2) = i+nun;
          Val_B(6*(k-1)+2) = 1.0;
          
          % A matrix entries
          
          I_A(4*(k-1)+1) = 2*(k-1)+1;
          J_A(4*(k-1)+1) = j1;
          Val_A(4*(k-1)+1) = 1-ColU(i+nc,j);

          I_A(4*(k-1)+2) = 2*(k-1)+1;
          J_A(4*(k-1)+2) = j2;
          Val_A(4*(k-1)+2) = ColU(i+nc,j);

          I_A(4*(k-1)+3) = 2*(k-1)+2;
          J_A(4*(k-1)+3) = j1+nc;
          Val_A(4*(k-1)+3) = 1-ColU(i+nc,j);

          I_A(4*(k-1)+4) = 2*(k-1)+2;
          J_A(4*(k-1)+4) = j2+nc;
          Val_A(4*(k-1)+4) = ColU(i+nc,j);
          
      end
      
  end

  N = sparse(I_N,J_N,Val_N,num_collisions,2*num_collisions);
  B = sparse(I_B,J_B,Val_B,2*num_collisions,2*nun);
  A = sparse(I_A,J_A,Val_A,2*num_collisions,2*nc); 