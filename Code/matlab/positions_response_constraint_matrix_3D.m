function [N,B,A,B1,B2] = positions_response_constraint_matrix_3D(F,Col,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,...
    EdgeCol,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2,V_fine,F_fine)
  % RESPONSE_CONSTRAINT_MATRIX_3D
  % response_constraint_matrix_3D(F,Col,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,...
  % EdgeCol,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2,V_fine,F_fine)
  %
  % Implements constraint matrices N^{un} and N^{c} as decribed in my
  % report.
  %
  % Input:
  %   F  N x 3 list of vertex indices
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
  %   (V_fine,F_fine) is the fine (constrained) mesh.
  % Output:
  %   N   (number of collisions)x(3*(number of collisions)) 
  %   matrix with entries being the entries of the collision normals
  %   B   (3*(number of collisions))x(3*num.vert. coarse mesh) barycentric
  %   coordinates of the vertices on the coarse mesh
  %   A   (3*(number of collisions))x(3*num.vert. fine mesh) barycentric
  %   coordinates of the vertices on the fine mesh


  num_collisions = full(sum(sum(Col)));
  edge_num_collisions = full(sum(sum(EdgeCol)));
  
  % number of vertices in the fine (constrained) mesh
  nc = size(V_fine,1);
  % number of vertices in the coarse (unconstrained) mesh
  nun = size(Col,1)-nc;
  
  [ColI,ColJ,~] = find(Col);
  [EdgeColI,EdgeColJ,~] = find(EdgeCol);
  
  I_N = ones(1,3*num_collisions);
  J_N = ones(1,3*num_collisions);
  Val_N = zeros(1,3*num_collisions);
  
  I_B = ones(1,12*num_collisions);
  J_B = ones(1,12*num_collisions);
  Val_B = zeros(1,12*num_collisions);
  
  I_A = ones(1,9*num_collisions);
  J_A = ones(1,9*num_collisions);
  Val_A = zeros(1,9*num_collisions);
  
  I_B1 = ones(1,12*num_collisions);
  J_B1 = ones(1,12*num_collisions);
  Val_B1 = zeros(1,12*num_collisions);
  
  I_B2 = ones(1,12*num_collisions);
  J_B2 = ones(1,12*num_collisions);
  Val_B2 = zeros(1,12*num_collisions);
  
  for k=1:num_collisions
      
      i = ColI(k);
      j = ColJ(k);
      
      % case 1: vertex i belongs to fine mesh, face j to the coarse mesh
      if i<=nc
          
          disp('entered case 1: vertex-face ');
                    
          j1 = F(j,1)-nc;
          j2 = F(j,2)-nc;
          j3 = F(j,3)-nc;

          % N matrix entries

          I_N(3*(k-1)+1) = k;
          J_N(3*(k-1)+1) = 3*(k-1)+1;
          Val_N(3*(k-1)+1) = ColNX(i,j);

          I_N(3*(k-1)+2) = k;
          J_N(3*(k-1)+2) = 3*(k-1)+2;
          Val_N(3*(k-1)+2) = ColNY(i,j);

          I_N(3*(k-1)+3) = k;
          J_N(3*(k-1)+3) = 3*(k-1)+3;
          Val_N(3*(k-1)+3) = ColNZ(i,j);

          % B matrix entries

          I_B(12*(k-1)+1) = 3*(k-1)+1;
          J_B(12*(k-1)+1) = j1;
          Val_B(12*(k-1)+1) = ColU1(i,j);

          I_B(12*(k-1)+2) = 3*(k-1)+1;
          J_B(12*(k-1)+2) = j2;
          Val_B(12*(k-1)+2) = ColU2(i,j);

          I_B(12*(k-1)+3) = 3*(k-1)+1;
          J_B(12*(k-1)+3) = j3;
          Val_B(12*(k-1)+3) = ColU3(i,j);

          I_B(12*(k-1)+4) = 3*(k-1)+2;
          J_B(12*(k-1)+4) = j1+nun;
          Val_B(12*(k-1)+4) = ColU1(i,j);

          I_B(12*(k-1)+5) = 3*(k-1)+2;
          J_B(12*(k-1)+5) = j2+nun;
          Val_B(12*(k-1)+5) = ColU2(i,j);

          I_B(12*(k-1)+6) = 3*(k-1)+2;
          J_B(12*(k-1)+6) = j3+nun;
          Val_B(12*(k-1)+6) = ColU3(i,j);

          I_B(12*(k-1)+7) = 3*(k-1)+3;
          J_B(12*(k-1)+7) = j1+2*nun;
          Val_B(12*(k-1)+7) = ColU1(i,j);

          I_B(12*(k-1)+8) = 3*(k-1)+3;
          J_B(12*(k-1)+8) = j2+2*nun;
          Val_B(12*(k-1)+8) = ColU2(i,j);

          I_B(12*(k-1)+9) = 3*(k-1)+3;
          J_B(12*(k-1)+9) = j3+2*nun;
          Val_B(12*(k-1)+9) = ColU3(i,j);

          % A matrix entries

          I_A(9*(k-1)+1) = 3*(k-1)+1;
          J_A(9*(k-1)+1) = i;
          Val_A(9*(k-1)+1) = 1.0;

          I_A(9*(k-1)+2) = 3*(k-1)+2;
          J_A(9*(k-1)+2) = i+nc;
          Val_A(9*(k-1)+2) = 1.0;

          I_A(9*(k-1)+3) = 3*(k-1)+3;
          J_A(9*(k-1)+3) = i+2*nc;
          Val_A(9*(k-1)+3) = 1.0;
                    
      % case 2: both vertex i and face j belong to coarse mesh    
      elseif (i>nc && j>size(F_fine,1))
          
          disp('entered case 2: vertex-face ');
          
%             ColU1 = ones(size(ColU1));
%             ColU2 = ones(size(ColU2));
%             ColU3 = ones(size(ColU3));
          
%           disp('vertex ');
%           i
%           disp('vertex in face');
%           F(j,1)
                    
          i = i-nc;
          j1 = F(j,1)-nc; j2 = F(j,2)-nc; j3 = F(j,3)-nc;
          
          % N matrix entries
          
          I_N(3*(k-1)+1) = k;
          J_N(3*(k-1)+1) = 3*(k-1)+1;
          Val_N(3*(k-1)+1) = ColNX(i+nc,j);
          
          I_N(3*(k-1)+2) = k;
          J_N(3*(k-1)+2) = 3*(k-1)+2;
          Val_N(3*(k-1)+2) = ColNY(i+nc,j);
          
          I_N(3*(k-1)+3) = k;
          J_N(3*(k-1)+3) = 3*(k-1)+3;
          Val_N(3*(k-1)+3) = ColNZ(i+nc,j);
          
          % B matrix entries

          I_B(12*(k-1)+1) = 3*(k-1)+1;
          J_B(12*(k-1)+1) = i;
          Val_B(12*(k-1)+1) = 1.0;

          I_B(12*(k-1)+2) = 3*(k-1)+2;
          J_B(12*(k-1)+2) = i+nun;
          Val_B(12*(k-1)+2) = 1.0;

          I_B(12*(k-1)+3) = 3*(k-1)+3;
          J_B(12*(k-1)+3) = i+2*nun;
          Val_B(12*(k-1)+3) = 1.0;
          
          I_B(12*(k-1)+4) = 3*(k-1)+1;
          J_B(12*(k-1)+4) = j1;
          Val_B(12*(k-1)+4) = -ColU1(i+nc,j);

          I_B(12*(k-1)+5) = 3*(k-1)+1;
          J_B(12*(k-1)+5) = j2;
          Val_B(12*(k-1)+5) = -ColU2(i+nc,j);

          I_B(12*(k-1)+6) = 3*(k-1)+1;
          J_B(12*(k-1)+6) = j3;
          Val_B(12*(k-1)+6) = -ColU3(i+nc,j);

          I_B(12*(k-1)+7) = 3*(k-1)+2;
          J_B(12*(k-1)+7) = j1+nun;
          Val_B(12*(k-1)+7) = -ColU1(i+nc,j);

          I_B(12*(k-1)+8) = 3*(k-1)+2;
          J_B(12*(k-1)+8) = j2+nun;
          Val_B(12*(k-1)+8) = -ColU2(i+nc,j);

          I_B(12*(k-1)+9) = 3*(k-1)+2;
          J_B(12*(k-1)+9) = j3+nun;
          Val_B(12*(k-1)+9) = -ColU3(i+nc,j);

          I_B(12*(k-1)+10) = 3*(k-1)+3;
          J_B(12*(k-1)+10) = j1+2*nun;
          Val_B(12*(k-1)+10) = -ColU1(i+nc,j);

          I_B(12*(k-1)+11) = 3*(k-1)+3;
          J_B(12*(k-1)+11) = j2+2*nun;
          Val_B(12*(k-1)+11) = -ColU2(i+nc,j);

          I_B(12*(k-1)+12) = 3*(k-1)+3;
          J_B(12*(k-1)+12) = j3+2*nun;
          Val_B(12*(k-1)+12) = -ColU3(i+nc,j);
          
          % debugging
          
          I_B1(12*(k-1)+1) = 3*(k-1)+1;
          J_B1(12*(k-1)+1) = i;
          Val_B1(12*(k-1)+1) = 1.0;

          I_B1(12*(k-1)+2) = 3*(k-1)+2;
          J_B1(12*(k-1)+2) = i+nun;
          Val_B1(12*(k-1)+2) = 1.0;

          I_B1(12*(k-1)+3) = 3*(k-1)+3;
          J_B1(12*(k-1)+3) = i+2*nun;
          Val_B1(12*(k-1)+3) = 1.0;
          
          I_B2(12*(k-1)+4) = 3*(k-1)+1;
          J_B2(12*(k-1)+4) = j1;
          Val_B2(12*(k-1)+4) = -ColU1(i+nc,j);

          I_B2(12*(k-1)+5) = 3*(k-1)+1;
          J_B2(12*(k-1)+5) = j2;
          Val_B2(12*(k-1)+5) = -ColU2(i+nc,j);

          I_B2(12*(k-1)+6) = 3*(k-1)+1;
          J_B2(12*(k-1)+6) = j3;
          Val_B2(12*(k-1)+6) = -ColU3(i+nc,j);

          I_B2(12*(k-1)+7) = 3*(k-1)+2;
          J_B2(12*(k-1)+7) = j1+nun;
          Val_B2(12*(k-1)+7) = -ColU1(i+nc,j);

          I_B2(12*(k-1)+8) = 3*(k-1)+2;
          J_B2(12*(k-1)+8) = j2+nun;
          Val_B2(12*(k-1)+8) = -ColU2(i+nc,j);

          I_B2(12*(k-1)+9) = 3*(k-1)+2;
          J_B2(12*(k-1)+9) = j3+nun;
          Val_B2(12*(k-1)+9) = -ColU3(i+nc,j);

          I_B2(12*(k-1)+10) = 3*(k-1)+3;
          J_B2(12*(k-1)+10) = j1+2*nun;
          Val_B2(12*(k-1)+10) = -ColU1(i+nc,j);

          I_B2(12*(k-1)+11) = 3*(k-1)+3;
          J_B2(12*(k-1)+11) = j2+2*nun;
          Val_B2(12*(k-1)+11) = -ColU2(i+nc,j);

          I_B2(12*(k-1)+12) = 3*(k-1)+3;
          J_B2(12*(k-1)+12) = j3+2*nun;
          Val_B2(12*(k-1)+12) = -ColU3(i+nc,j);
          
%           xb = ColU1(i,j)*V0(F(j,1),:)+ColU2(i,j)*V0(F(j,2),:)+ColU3(i,j)*V0(F(j,3),:);
%           xa = V0(i,:);
%           
%           dot(,[ColNX(i+nc,j) ColNY(i+nc,j) ColNZ(i+nc,j)])
      
      % case 3: face j belongs to fine mesh, vertex i belongs to coarse
      else
          
          disp('entered case 3: vertex-face ');
                              
          i = i-nc;
          j1 = F(j,1); j2 = F(j,2); j3 = F(j,3);
          
          % N matrix entries

          I_N(3*(k-1)+1) = k;
          J_N(3*(k-1)+1) = 3*(k-1)+1;
          Val_N(3*(k-1)+1) = ColNX(i+nc,j);

          I_N(3*(k-1)+2) = k;
          J_N(3*(k-1)+2) = 3*(k-1)+2;
          Val_N(3*(k-1)+2) = ColNY(i+nc,j);

          I_N(3*(k-1)+3) = k;
          J_N(3*(k-1)+3) = 3*(k-1)+3;
          Val_N(3*(k-1)+3) = ColNZ(i+nc,j);
          
          % B matrix entries

          I_B(12*(k-1)+1) = 3*(k-1)+1;
          J_B(12*(k-1)+1) = i;
          Val_B(12*(k-1)+1) = 1.0;

          I_B(12*(k-1)+2) = 3*(k-1)+2;
          J_B(12*(k-1)+2) = i+nun;
          Val_B(12*(k-1)+2) = 1.0;

          I_B(12*(k-1)+3) = 3*(k-1)+3;
          J_B(12*(k-1)+3) = i+2*nun;
          Val_B(12*(k-1)+3) = 1.0;
          
          % A matrix entries
          
          I_A(9*(k-1)+1) = 3*(k-1)+1;
          J_A(9*(k-1)+1) = j1;
          Val_A(9*(k-1)+1) = ColU1(i+nc,j);

          I_A(9*(k-1)+2) = 3*(k-1)+1;
          J_A(9*(k-1)+2) = j2;
          Val_A(9*(k-1)+2) = ColU2(i+nc,j);

          I_A(9*(k-1)+3) = 3*(k-1)+1;
          J_A(9*(k-1)+3) = j3;
          Val_A(9*(k-1)+3) = ColU3(i+nc,j);

          I_A(9*(k-1)+4) = 3*(k-1)+2;
          J_A(9*(k-1)+4) = j1+nc;
          Val_A(9*(k-1)+4) = ColU1(i+nc,j);

          I_A(9*(k-1)+5) = 3*(k-1)+2;
          J_A(9*(k-1)+5) = j2+nc;
          Val_A(9*(k-1)+5) = ColU2(i+nc,j);

          I_A(9*(k-1)+6) = 3*(k-1)+2;
          J_A(9*(k-1)+6) = j3+nc;
          Val_A(9*(k-1)+6) = ColU3(i+nc,j);

          I_A(9*(k-1)+7) = 3*(k-1)+3;
          J_A(9*(k-1)+7) = j1+2*nc;
          Val_A(9*(k-1)+7) = ColU1(i+nc,j);

          I_A(9*(k-1)+8) = 3*(k-1)+3;
          J_A(9*(k-1)+8) = j2+2*nc;
          Val_A(9*(k-1)+8) = ColU2(i+nc,j);

          I_A(9*(k-1)+9) = 3*(k-1)+3;
          J_A(9*(k-1)+9) = j3+2*nc;
          Val_A(9*(k-1)+9) = ColU3(i+nc,j);
          
      end
      
  end
  
  N = sparse(I_N,J_N,Val_N,num_collisions,3*num_collisions+3*edge_num_collisions);
  B = sparse(I_B,J_B,Val_B,3*num_collisions,3*nun);
  A = sparse(I_A,J_A,Val_A,3*num_collisions,3*nc);
  B1 = sparse(I_B1,J_B1,Val_B1,3*num_collisions,3*nun);
  B2 = sparse(I_B2,J_B2,Val_B2,3*num_collisions,3*nun);
%   
%   Val_B2
  
  % now edge-edge collisions
  
  % compute all edges
  EE_all = edges(F);
  
  I_N = ones(1,3*edge_num_collisions);
  J_N = ones(1,3*edge_num_collisions);
  Val_N = zeros(1,3*edge_num_collisions);
  
  I_B = ones(1,12*edge_num_collisions);
  J_B = ones(1,12*edge_num_collisions);
  Val_B = zeros(1,12*edge_num_collisions);
  
  I_A = ones(1,12*edge_num_collisions);
  J_A = ones(1,12*edge_num_collisions);
  Val_A = zeros(1,12*edge_num_collisions);
  
  for k=1:edge_num_collisions
      
      i = EdgeColI(k);
      j = EdgeColJ(k);
      u1 = EdgeColU1(i,j); u2 = EdgeColU2(i,j);
      
      i1 = EE_all(i,1);
      i2 = EE_all(i,2);
      j1 = EE_all(j,1);
      j2 = EE_all(j,2);
      
      % Case 1: edge j belongs to fine mesh and edge i belongs to coarse
      if (j1<=size(V_fine,1))
          
          disp('entered case 1: edge-edge ');
                    
          i1 = EE_all(i,1)-nc; i2 = EE_all(i,2)-nc; 
          j1 = EE_all(j,1); j2 = EE_all(j,2); 
      
          % N matrix entries

          I_N(3*(k-1)+1) = k;
          J_N(3*(k-1)+1) = 3*num_collisions+3*(k-1)+1;
          Val_N(3*(k-1)+1) = EdgeColNX(i,j);

          I_N(3*(k-1)+2) = k;
          J_N(3*(k-1)+2) = 3*num_collisions+3*(k-1)+2;
          Val_N(3*(k-1)+2) = EdgeColNY(i,j);

          I_N(3*(k-1)+3) = k;
          J_N(3*(k-1)+3) = 3*num_collisions+3*(k-1)+3;
          Val_N(3*(k-1)+3) = EdgeColNZ(i,j);

          % B matrix entries

          I_B(12*(k-1)+1) = 3*(k-1)+1;
          J_B(12*(k-1)+1) = i1;
          Val_B(12*(k-1)+1) = 1-u1;

          I_B(12*(k-1)+2) = 3*(k-1)+1;
          J_B(12*(k-1)+2) = i2;
          Val_B(12*(k-1)+2) = u1;

          I_B(12*(k-1)+3) = 3*(k-1)+2;
          J_B(12*(k-1)+3) = i1+nun;
          Val_B(12*(k-1)+3) = 1-u1;

          I_B(12*(k-1)+4) = 3*(k-1)+2;
          J_B(12*(k-1)+4) = i2+nun;
          Val_B(12*(k-1)+4) = u1;

          I_B(12*(k-1)+5) = 3*(k-1)+3;
          J_B(12*(k-1)+5) = i1+2*nun;
          Val_B(12*(k-1)+5) = 1-u1;

          I_B(12*(k-1)+6) = 3*(k-1)+3;
          J_B(12*(k-1)+6) = i2+2*nun;
          Val_B(12*(k-1)+6) = u1;

          % A matrix entries

          I_A(12*(k-1)+1) = 3*(k-1)+1;
          J_A(12*(k-1)+1) = j1;
          Val_A(12*(k-1)+1) = 1-u2;

          I_A(12*(k-1)+2) = 3*(k-1)+1;
          J_A(12*(k-1)+2) = j2;
          Val_A(12*(k-1)+2) = u2;

          I_A(12*(k-1)+3) = 3*(k-1)+2;
          J_A(12*(k-1)+3) = j1+nc;
          Val_A(12*(k-1)+3) = 1-u2;

          I_A(12*(k-1)+4) = 3*(k-1)+2;
          J_A(12*(k-1)+4) = j2+nc;
          Val_A(12*(k-1)+4) = u2;

          I_A(12*(k-1)+5) = 3*(k-1)+3;
          J_A(12*(k-1)+5) = j1+2*nc;
          Val_A(12*(k-1)+5) = 1-u2;

          I_A(12*(k-1)+6) = 3*(k-1)+3;
          J_A(12*(k-1)+6) = j2+2*nc;
          Val_A(12*(k-1)+6) = u2;
          
      % Case 2: both edge j and edge i belongs to coarse mesh   
      elseif (j1>size(V_fine,1) && i1>size(V_fine,1))
          
          disp('entered case 2: edge-edge ');
          
%           disp('vertex in edge 1 ')
%           EE_all(i,1)
%           
%           disp('vertex in edge 2 ')
%           EE_all(j,1)
          
          i1 = EE_all(i,1)-nc; i2 = EE_all(i,2)-nc; 
          j1 = EE_all(j,1)-nc; j2 = EE_all(j,2)-nc; 
      
          % N matrix entries

          I_N(3*(k-1)+1) = k;
          J_N(3*(k-1)+1) = 3*num_collisions+3*(k-1)+1;
          Val_N(3*(k-1)+1) = EdgeColNX(i,j);

          I_N(3*(k-1)+2) = k;
          J_N(3*(k-1)+2) = 3*num_collisions+3*(k-1)+2;
          Val_N(3*(k-1)+2) = EdgeColNY(i,j);

          I_N(3*(k-1)+3) = k;
          J_N(3*(k-1)+3) = 3*num_collisions+3*(k-1)+3;
          Val_N(3*(k-1)+3) = EdgeColNZ(i,j);
          
          % B matrix entries

          I_B(12*(k-1)+1) = 3*(k-1)+1;
          J_B(12*(k-1)+1) = j1;
          Val_B(12*(k-1)+1) = 1-u2;

          I_B(12*(k-1)+2) = 3*(k-1)+1;
          J_B(12*(k-1)+2) = j2;
          Val_B(12*(k-1)+2) = u2;

          I_B(12*(k-1)+3) = 3*(k-1)+2;
          J_B(12*(k-1)+3) = j1+nun;
          Val_B(12*(k-1)+3) = 1-u2;

          I_B(12*(k-1)+4) = 3*(k-1)+2;
          J_B(12*(k-1)+4) = j2+nun;
          Val_B(12*(k-1)+4) = u2;

          I_B(12*(k-1)+5) = 3*(k-1)+3;
          J_B(12*(k-1)+5) = j1+2*nun;
          Val_B(12*(k-1)+5) = 1-u2;

          I_B(12*(k-1)+6) = 3*(k-1)+3;
          J_B(12*(k-1)+6) = j2+2*nun;
          Val_B(12*(k-1)+6) = u2;
          
          I_B(12*(k-1)+7) = 3*(k-1)+1;
          J_B(12*(k-1)+7) = i1;
          Val_B(12*(k-1)+7) = -(1-u1);

          I_B(12*(k-1)+8) = 3*(k-1)+1;
          J_B(12*(k-1)+8) = i2;
          Val_B(12*(k-1)+8) = -u1;

          I_B(12*(k-1)+9) = 3*(k-1)+2;
          J_B(12*(k-1)+9) = i1+nun;
          Val_B(12*(k-1)+9) = -(1-u1);

          I_B(12*(k-1)+10) = 3*(k-1)+2;
          J_B(12*(k-1)+10) = i2+nun;
          Val_B(12*(k-1)+10) = -u1;

          I_B(12*(k-1)+11) = 3*(k-1)+3;
          J_B(12*(k-1)+11) = i1+2*nun;
          Val_B(12*(k-1)+11) = -(1-u1);

          I_B(12*(k-1)+12) = 3*(k-1)+3;
          J_B(12*(k-1)+12) = i2+2*nun;
          Val_B(12*(k-1)+12) = -u1;
          
      
      % Case 3: edge j belongs to coarse mesh and edge i belongs to fine
      else
          
          disp('entered case 3: edge-edge ');
          
          i1 = EE_all(i,1); i2 = EE_all(i,2); 
          j1 = EE_all(j,1)-nc; j2 = EE_all(j,2)-nc; 
      
          % N matrix entries

          I_N(3*(k-1)+1) = k;
          J_N(3*(k-1)+1) = 3*num_collisions+3*(k-1)+1;
          Val_N(3*(k-1)+1) = EdgeColNX(i,j);

          I_N(3*(k-1)+2) = k;
          J_N(3*(k-1)+2) = 3*num_collisions+3*(k-1)+2;
          Val_N(3*(k-1)+2) = EdgeColNY(i,j);

          I_N(3*(k-1)+3) = k;
          J_N(3*(k-1)+3) = 3*num_collisions+3*(k-1)+3;
          Val_N(3*(k-1)+3) = EdgeColNZ(i,j);
          
          % B matrix entries

          I_B(12*(k-1)+1) = 3*(k-1)+1;
          J_B(12*(k-1)+1) = j1;
          Val_B(12*(k-1)+1) = 1-u2;

          I_B(12*(k-1)+2) = 3*(k-1)+1;
          J_B(12*(k-1)+2) = j2;
          Val_B(12*(k-1)+2) = u2;

          I_B(12*(k-1)+3) = 3*(k-1)+2;
          J_B(12*(k-1)+3) = j1+nun;
          Val_B(12*(k-1)+3) = 1-u2;

          I_B(12*(k-1)+4) = 3*(k-1)+2;
          J_B(12*(k-1)+4) = j2+nun;
          Val_B(12*(k-1)+4) = u2;

          I_B(12*(k-1)+5) = 3*(k-1)+3;
          J_B(12*(k-1)+5) = j1+2*nun;
          Val_B(12*(k-1)+5) = 1-u2;

          I_B(12*(k-1)+6) = 3*(k-1)+3;
          J_B(12*(k-1)+6) = j2+2*nun;
          Val_B(12*(k-1)+6) = u2;
          
          % A matrix entries

          I_A(12*(k-1)+1) = 3*(k-1)+1;
          J_A(12*(k-1)+1) = i1;
          Val_A(12*(k-1)+1) = 1-u1;

          I_A(12*(k-1)+2) = 3*(k-1)+1;
          J_A(12*(k-1)+2) = i2;
          Val_A(12*(k-1)+2) = u1;

          I_A(12*(k-1)+3) = 3*(k-1)+2;
          J_A(12*(k-1)+3) = i1+nc;
          Val_A(12*(k-1)+3) = 1-u1;

          I_A(12*(k-1)+4) = 3*(k-1)+2;
          J_A(12*(k-1)+4) = i2+nc;
          Val_A(12*(k-1)+4) = u1;

          I_A(12*(k-1)+5) = 3*(k-1)+3;
          J_A(12*(k-1)+5) = i1+2*nc;
          Val_A(12*(k-1)+5) = 1-u1;

          I_A(12*(k-1)+6) = 3*(k-1)+3;
          J_A(12*(k-1)+6) = i2+2*nc;
          Val_A(12*(k-1)+6) = u1;
          
      end
      
      
  end
  
  if (edge_num_collisions>0)
      
      N = [N; sparse(I_N,J_N,Val_N,edge_num_collisions,3*num_collisions+3*edge_num_collisions)];
      B = [B; sparse(I_B,J_B,Val_B,3*edge_num_collisions,3*nun)];
      A = [A; sparse(I_A,J_A,Val_A,3*edge_num_collisions,3*nc)];
  end