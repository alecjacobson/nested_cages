function [Col,ColX,ColY] = collision_detection_2d(V0,V1,E)
  % COLLISION_DETECTION_2D
  % collision_detection_2d(V0,V1,E,varagin)
  %
  % Implements vertex-edge collision detection as described in "Collision
  % and Self-Collision Handling..." by Provot 97 (section 1.1, changing
  % triangle by edge)
  %
  % Input:
  %   V0  M x 2 list of polygon vertex positions at the initial time t0
  %   V1  M x 2 list of polygon vertex positions at the next time
  %   t0+delta_t
  %   E0  2 x N list of vertex indices that form the edges for both meshes 
  % Output:
  %   Col   MxN matrix such that Col(i,j) == 1 iff there is a collision between
  %     vertex i and edges j
  %   ColX   MxN matrix such that ColY(i,j) is the y-coordinates
  %     2D position where vertex i collides with edge j. If there is no
  %   ColY   MxN matrix such that ColY(i,j) is the y-coordinates
  %     2D position where vertex i collides with edge j. If there is no
  
  % rotated edges at time t0
  Edges_0 = V0(E(1,:),:)-V0(E(2,:),:);
  Normals_0 = [-Edges_0(:,2) Edges_0(:,1)];
  
  % rotated edges at time t0+delta_t
  Edges_1 = V1(E(1,:),:)-V1(E(2,:),:);
  Normals_1 = [-Edges_1(:,2) Edges_1(:,1)];
  
  % rotated edges at time t0
  Edges_0 = V0(E(1,:),:)-V0(E(2,:),:);
  Normals_0 = [-Edges_0(:,2) Edges_0(:,1)];
  
  % rotated edges at time t0+delta_t
  Edges_1 = V1(E(1,:),:)-V1(E(2,:),:);
  Normals_1 = [-Edges_1(:,2) Edges_1(:,1)];

  % Initialize output collisions as +Infinity
  Col = sparse(size(V0,1),size(E,2));
  ColX = sparse(size(V0,1),size(E,2));
  ColY = sparse(size(V0,1),size(E,2));
  % Initialize coefficients of the polynomial to be solved as zeros
  p = zeros(1,3);
  
  % loop over all edges j and vertices i
  for j=1:size(E,2)
      
      for i=setdiff(1:size(V0,1),[E(1,j) E(2,j)])
          
          % precompute inner products that will be needed
          inner00 = dot(Normals_0(j,:),V0(i,:)-V0(E(1,j),:));
          inner01 = dot(Normals_0(j,:),V1(i,:)-V1(E(1,j),:));
          inner10 = dot(Normals_1(j,:),V0(i,:)-V0(E(1,j),:));
          inner11 = dot(Normals_1(j,:),V1(i,:)-V1(E(1,j),:));
          
          % calculate coefficients of the polynomial
          p(1) = inner00-inner01-inner10+inner11;
          p(2) = -2*inner00+inner01+inner10;
          p(3) = inner00;
          
          % find its roots
          sols = roots(p);
          
          % select roots that belong to [0,1]. If more than one,
          % select the smallest one
          t = sols(find((sols>=0).*(sols<=1)));
          t = min(t);
          
          % if returned some t from previous selection
          if(size(t,1)==1)
              
              % compute all necessary points
              P_t = (1-t)*V0(i,:)+t*V1(i,:);
              A_t = (1-t)*V0(E(1,j),:)+t*V1(E(1,j),:);
              B_t = (1-t)*V0(E(2,j),:)+t*V1(E(2,j),:);

              % determine if P lies on AB edge (it does if u below belongs
              % to [0,1])
              [~,idx] = max(abs(B_t-A_t));
              u = (P_t(idx)-A_t(idx))/(B_t(idx)-A_t(idx));

              % if it u belongs to [0,1], output colliding point
              if (u>=0&&u<=1)
                       Col(i, j) = 1;
                       ColX(i,j) = P_t(1);
                       ColY(i,j) = P_t(2);
              end
              
          end
          
      end
      
  end
