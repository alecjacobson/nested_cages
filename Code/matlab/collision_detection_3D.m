function [Col,ColX,ColY,ColZ,ColNX,ColNY,ColNZ,ColU1,ColU2,ColU3,ColT,...
    EdgeCol,EdgeColX,EdgeColY,EdgeColZ,EdgeColNX,EdgeColNY,EdgeColNZ,EdgeColU1,EdgeColU2,EdgeColT,...
    ColXA,ColYA,ColZA,ColXB,ColYB,ColZB]...
    = collision_detection_3D(V0,V1,F,V_fine,F_fine,varargin)
  % COLLISION_DETECTION_3D
  % collision_detection_3d(V0,V1,F,varargin)
  %
  % Vertex-face and edge-edge collision detection
  %
  % Input:
  %   V0  M x 3 list of mesh vertex positions at the initial time t0
  %   V1  M x 3 list of mesh vertex positions at the next time
  %   t0+delta_t
  %   F0  Nx3 list of vertex indices that form the faces for both meshes
  %   (V_fine,F_fine): fine mesh
  %   Optional:
  %     'last_only' test only if last point of the mesh collides
  % Output:
  %   Col   MxN matrix such that Col(i,j) == 1 iff there is a collision between
  %     vertex i and face j
  %   ColX   MxN matrix such that ColX(i,j) is the x-coordinates
  %     3D position where vertex i collides with face j. 
  %   ColY   MxN matrix such that ColY(i,j) is the y-coordinates
  %     3D position where vertex i collides with face j. 
  %   ColZ   MxN matrix such that ColZ(i,j) is the z-coordinates
  %     3D position where vertex i collides with face j. 
  %   ColNX   MxN matrix such that ColNX(i,j) is the x-coordinate
  %     of the collision normal of vertex i collides with face j. 
  %   ColNY   MxN matrix such that ColNY(i,j) is the y-coordinate
  %     of the collision normal of vertex i collides with face j. 
  %   ColNZ   MxN matrix such that ColNZ(i,j) is the z-coordinate
  %     of the collision normal of vertex i collides with face j.
  %   ColU1   MxN matrix such that ColU(i,j) is the barycentric 
  %   coordinate of the position where vertex i collides with face j,
  %   w.r.t. the 1st vertex of face j. 
  %   ColU2   MxN matrix such that ColU(i,j) is the barycentric 
  %   coordinate of the position where vertex i collides with face j,
  %   w.r.t. the 2nd vertex of face j.
  %   ColU3   MxN matrix such that ColU(i,j) is the barycentric 
  %   coordinate of the position where vertex i collides with face j,
  %   w.r.t. the 3rd vertex of face j.
  %   ColT   MxN matrix such that ColT(i,j) is the time instant when
  %   vertex i collides with edge j.
  %   EdgeCol   (#edges)x(#edges) matrix such that Col(i,j) == 1 
  %     iff there is a collision between edge i and edge j
  %   EdgeColX   (#edges)x(#edges) matrix such that ColX(i,j) 
  %   is the x-coordinate of the position position where 
  %   edge i collides with edge j. 
  %   EdgeColY   (#edges)x(#edges) matrix such that ColX(i,j) 
  %   is the x-coordinate of the position position where 
  %   edge i collides with edge j. 
  %   EdgeColZ   (#edges)x(#edges) matrix such that ColX(i,j) 
  %   is the x-coordinate of the position position where 
  %   edge i collides with edge j. 
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
  
  last_only = 0;
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'last_only'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              last_only = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % calculate normals based on laplacian of the positions
  L = cotmatrix(V1,F);
  M = massmatrix(V1,F,'voronoi');
  Lap_Embedding = M*L*V1;
  Lap_Normals = -normalizerow(Lap_Embedding);
  
  % Initialize output collisions as zero
  Col = sparse(size(V0,1),size(F,1));
  ColX = sparse(size(V0,1),size(F,1));
  ColY = sparse(size(V0,1),size(F,1));
  ColZ = sparse(size(V0,1),size(F,1));
  ColNX = sparse(size(V0,1),size(F,1));
  ColNY = sparse(size(V0,1),size(F,1));
  ColNZ = sparse(size(V0,1),size(F,1));
  ColU1 = sparse(size(V0,1),size(F,1));
  ColU2 = sparse(size(V0,1),size(F,1));
  ColU3 = sparse(size(V0,1),size(F,1));
  ColT = sparse(size(V0,1),size(F,1));
  
  ColXA = sparse(size(V0,1),size(F,1));
  ColYA = sparse(size(V0,1),size(F,1));
  ColZA = sparse(size(V0,1),size(F,1));
  ColXB = sparse(size(V0,1),size(F,1));
  ColYB = sparse(size(V0,1),size(F,1));
  ColZB = sparse(size(V0,1),size(F,1));
    
  % collision detection
  [VF,EE] = ccd_volume_mex(V0,V1,F,size(V_fine,1));
  
  % loop over collisions to collect the data we need for the simulation
  for k=1:size(VF,1)
      
      % find vertex and face
      i = VF(k,1);

      [~,j]=ismember([VF(k,2) VF(k,3) VF(k,4)],F,'rows');
      
      Col(i,j) = 1;
      
      % time instant
      ColT(i,j) = VF(k,5);
      
      % positions
      ColX(i,j) = (1-VF(k,5))*V0(i,1)+(VF(k,5))*V1(i,1);
      ColY(i,j) = (1-VF(k,5))*V0(i,2)+(VF(k,5))*V1(i,2);
      ColZ(i,j) = (1-VF(k,5))*V0(i,3)+(VF(k,5))*V1(i,3);

      if (i<size(V_fine,1))
%           % if vertex belogns to fine mesh, take Laplacian Normals
%           N1 = Lap_Normals(i,:);

%           % fine mesh normals
%           idxs = find(F==i);
%           N1 = zeros(1,3);
%           for b=1:size(idxs)
%               f0(b) = mod(idxs(b),size(F,1));
%               if f0(b)==0
%                   f0(b) = size(F,1);
%               end
%               N1 = N1 + normals([V1(F(f0(b),1),:);V1(F(f0(b),2),:);V1(F(f0(b),3),:)],[1 2 3]);
%           end

          % coarse mesh normals
          N1 = normals([V0(F(j,1),:);V0(F(j,2),:);V0(F(j,3),:)],[1 2 3]);
          
          N1 = N1/norm(N1);
      else
          % if face belong to fine mesh, take face normals
          N1 = normals([V0(F(j,1),:);V0(F(j,2),:);V0(F(j,3),:)],[1 2 3]);      
          N1 = N1/norm(N1);
      end
      
      P = [ColX(i,j) ColY(i,j) ColZ(i,j)];
            
      % barycentric coordinates
      area1 = doublearea([(1-VF(k,5))*V0(F(j,2),:)+(VF(k,5))*V1(F(j,2),:);(1-VF(k,5))*V0(F(j,3),:)+(VF(k,5))*V1(F(j,3),:);P],[1 2 3])/2;
      area2 = doublearea([(1-VF(k,5))*V0(F(j,3),:)+(VF(k,5))*V1(F(j,3),:);(1-VF(k,5))*V0(F(j,1),:)+(VF(k,5))*V1(F(j,1),:);P],[1 2 3])/2;
      area3 = doublearea([(1-VF(k,5))*V0(F(j,1),:)+(VF(k,5))*V1(F(j,1),:);(1-VF(k,5))*V0(F(j,2),:)+(VF(k,5))*V1(F(j,2),:);P],[1 2 3])/2;
      sum_areas = area1+area2+area3;
      ColU1(i,j) = area1/sum_areas;
      ColU2(i,j) = area2/sum_areas;
      ColU3(i,j) = area3/sum_areas;
      
      
%       disp('position from collision code')
%       P

%       disp('barycentric coordinates from collision code')
%       ColU1(i,j)
%       ColU2(i,j)
%       ColU3(i,j)
%       disp('positions from collision code')
%       V0(F(j,1),:)
%       V0(F(j,2),:)
%       V0(F(j,3),:)
%       disp('face vertex indices from collision code')
%       F(j,1)-size(V_fine,1)
%       F(j,2)-size(V_fine,1)
%       F(j,3)-size(V_fine,1)
      
%       disp('vertex index')
%       i-size(V_fine,1)
%       
%       disp('normal from inside collision')
%       N1
      
      xb = V0(i,:);
      xa = ColU1(i,j)*V0(F(j,1),:)+ColU2(i,j)*V0(F(j,2),:)+ColU3(i,j)*V0(F(j,3),:);
      
%       % barycentric coordinates
%       disp('xb from collision code')
%       xb
%       disp('xa from collision code')
%       xa
%       disp('xb-xa from collision code')
%       xb-xa
      
%       disp('xb')
%       xb
%       disp('V0(n_fine+131,:) from inside the collision code')
%       V0(size(V_fine,1)+131,:)
      
%       disp('vertex at end of time')
%       V1(i,:)
      
%       disp('dot product at time 0');
%       dot(xb-xa,N1)
%       disp('distance at time 0 - vertex-face')
%       norm(P-V0(i,:))
      if (dot(xb-xa,N1)<0 && i>size(V_fine,1) && F(j,1)>size(V_fine,1))
          disp('WARNING: dot(xb-xa,N1)>0 and both vertices belong to coarse mesh')
          N1 = -N1;
      end
      xb = ColU1(i,j)*V1(F(j,1),:)+ColU2(i,j)*V1(F(j,2),:)+ColU3(i,j)*V1(F(j,3),:);
      xa = V1(i,:); 
%       disp('dot product at time 1');
%       dot(xb-xa,N1)
%       disp('vertex indices: vertex and vertex on face');
%       i
%       F(j,1)

      ColXA(i,j) = xa(1);
      ColYA(i,j) = xa(2);
      ColZA(i,j) = xa(3);
      
      ColXB(i,j) = xb(1);
      ColYB(i,j) = xb(2);
      ColZB(i,j) = xb(3);

      ColNX(i,j) = N1(1);
      ColNY(i,j) = N1(2);
      ColNZ(i,j) = N1(3);  
      
  end
  
%     %Initialize output collisions as zero
%   Col = sparse(size(V0,1),size(F,1));
%   ColX = sparse(size(V0,1),size(F,1));
%   ColY = sparse(size(V0,1),size(F,1));
%   ColZ = sparse(size(V0,1),size(F,1));
%   ColNX = sparse(size(V0,1),size(F,1));
%   ColNY = sparse(size(V0,1),size(F,1));
%   ColNZ = sparse(size(V0,1),size(F,1));
%   ColU1 = sparse(size(V0,1),size(F,1));
%   ColU2 = sparse(size(V0,1),size(F,1));
%   ColU3 = sparse(size(V0,1),size(F,1));
%   ColT = sparse(size(V0,1),size(F,1));
  
  % compute all edges
  EE_all = edges(F);
  
  % Initialize output collisions as zero
  EdgeCol = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColX = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColY = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColZ = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColNX = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColNY = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColNZ = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColU1 = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColU2 = sparse(size(EE_all,1),size(EE_all,1));
  EdgeColT = sparse(size(EE_all,1),size(EE_all,1));
  
  for k=1:size(EE,1)
      
    [~,i]=ismember([EE(k,2) EE(k,1)],EE_all,'rows');
    [~,j]=ismember([EE(k,4) EE(k,3)],EE_all,'rows');
    
    EdgeCol(i,j) = 1;
    t_col = EE(k,5);
    EdgeColT(i,j) = t_col;
    
    A1t = (1-t_col)*V0(EE_all(i,1),:)+t_col*V1(EE_all(i,1),:);
    B1t = (1-t_col)*V0(EE_all(i,2),:)+t_col*V1(EE_all(i,2),:);
    A2t = (1-t_col)*V0(EE_all(j,1),:)+t_col*V1(EE_all(j,1),:);
    B2t = (1-t_col)*V0(EE_all(j,2),:)+t_col*V1(EE_all(j,2),:);
    
    u1 = (A1t(2)*A2t(1) - A1t(3)*A2t(1) - A1t(1)*A2t(2) + A1t(3)*A2t(2) + A1t(1)*A2t(3) - ...
    A1t(2)*A2t(3) - A1t(2)*B2t(1) + A1t(3)*B2t(1) + A2t(2)*B2t(1) - A2t(3)*B2t(1) + ...
    A1t(1)*B2t(2) - A1t(3)*B2t(2) - A2t(1)*B2t(2) + A2t(3)*B2t(2) - A1t(1)*B2t(3) + ...
    A1t(2)*B2t(3) + A2t(1)*B2t(3) - A2t(2)*B2t(3))/(A1t(2)*A2t(1) - A1t(3)*A2t(1) - ...
    A1t(1)*A2t(2) + A1t(3)*A2t(2) + A1t(1)*A2t(3) - A1t(2)*A2t(3) + A2t(2)*B1t(1) - ...
    A2t(3)*B1t(1) - A2t(1)*B1t(2) + A2t(3)*B1t(2) + A2t(1)*B1t(3) - A2t(2)*B1t(3) - ...
    A1t(2)*B2t(1) + A1t(3)*B2t(1) + B1t(2)*B2t(1) - B1t(3)*B2t(1) + A1t(1)*B2t(2) - ...
    A1t(3)*B2t(2) - B1t(1)*B2t(2) + B1t(3)*B2t(2) - A1t(1)*B2t(3) + A1t(2)*B2t(3) + ...
    B1t(1)*B2t(3) - B1t(2)*B2t(3));

    EdgeColU1(i,j) = u1;
    
    P = (1-u1)*A1t+u1*B1t;
    EdgeColX(i,j) = P(1);
    EdgeColY(i,j) = P(2);
    EdgeColZ(i,j) = P(3);
    
    u2 = 1-(norm(P-B2t)/norm(A2t-B2t));
    EdgeColU2(i,j) = u2;
    
    % to calcullate the collision normals
    A10 = V0(EE_all(i,1),:);
    B10 = V0(EE_all(i,2),:);
    A20 = V0(EE_all(j,1),:);
    B20 = V0(EE_all(j,2),:);
    
    % to calcullate the collision normals
    A11 = V1(EE_all(i,1),:);
    B11 = V1(EE_all(i,2),:);
    A21 = V1(EE_all(j,1),:);
    B21 = V1(EE_all(j,2),:);
    
%     % fine mesh normals
%     if (i<j)
%         f = find(sum(ismember(F,[EE(k,2) EE(k,1)]),2)>=2);
%     else
%         f = find(sum(ismember(F,[EE(k,4) EE(k,3)]),2)>=2);
%     end
% 
%     N1 = normals([V1(F(f(1),1),:);V1(F(f(1),2),:);V1(F(f(1),3),:)],[1 2 3]);
%     N2 = normals([V1(F(f(2),1),:);V1(F(f(2),2),:);V1(F(f(2),3),:)],[1 2 3]);
%     N = (N1+N2)/norm(N1+N2);
    
      % coarse mesh normals
    if (i<j)
        f = find(sum(ismember(F,[EE(k,2) EE(k,1)]),2)>=2);
        disp('ENTERED IF');
        
%         disp('faces on the coarse mesh')
%         g = find(sum(ismember(F,[EE(k,2) EE(k,1)]),2)>=2)';
%         g - size(F_fine,1)
%         disp('one ring on the coarse mesh')
%         face_idx = [];
%         for z = 1:size(g,2)
%             face_idx = [face_idx (mod(find(F==F(g(z),1))-1,size(F,1))+1)' (mod(find(F==F(g(z),2))-1,size(F,1))+1)' (mod(find(F==F(g(z),3))-1,size(F,1))+1)'];
%         end
%         g = unique(face_idx);
%         g - size(F_fine,1)
%         disp('two ring on the coarse mesh')
%         face_idx = [];
%         for z = 1:size(g,2)
%             face_idx = [face_idx (mod(find(F==F(g(z),1))-1,size(F,1))+1)' (mod(find(F==F(g(z),2))-1,size(F,1))+1)' (mod(find(F==F(g(z),3))-1,size(F,1))+1)'];
%         end
%         g = unique(face_idx);
%         g - size(F_fine,1)
%         disp('three ring on the coarse mesh')
%         face_idx = [];
%         for z = 1:size(g,2)
%             face_idx = [face_idx (mod(find(F==F(g(z),1))-1,size(F,1))+1)' (mod(find(F==F(g(z),2))-1,size(F,1))+1)' (mod(find(F==F(g(z),3))-1,size(F,1))+1)'];
%         end
%         g = unique(face_idx);
%         g - size(F_fine,1)
%         
%         disp('faces on the fine mesh')
%         g = find(sum(ismember(F,[EE(k,4) EE(k,3)]),2)>=2)';
%         g
%         disp('one ring on the fine mesh')
%         face_idx = [];
%         for z = 1:size(g,2)
%             face_idx = [face_idx (mod(find(F==F(g(z),1))-1,size(F,1))+1)' (mod(find(F==F(g(z),2))-1,size(F,1))+1)' (mod(find(F==F(g(z),3))-1,size(F,1))+1)'];
%         end
%         g = unique(face_idx);
%         g
%         disp('two ring on the fine mesh')
%         face_idx = [];
%         for z = 1:size(g,2)
%             face_idx = [face_idx (mod(find(F==F(g(z),1))-1,size(F,1))+1)' (mod(find(F==F(g(z),2))-1,size(F,1))+1)' (mod(find(F==F(g(z),3))-1,size(F,1))+1)'];
%         end
%         g = unique(face_idx);
%         g
%         disp('three ring on the fine mesh')
%         face_idx = [];
%         for z = 1:size(g,2)
%             face_idx = [face_idx (mod(find(F==F(g(z),1))-1,size(F,1))+1)' (mod(find(F==F(g(z),2))-1,size(F,1))+1)' (mod(find(F==F(g(z),3))-1,size(F,1))+1)'];
%         end
%         g = unique(face_idx);
%         g
    else
        f = find(sum(ismember(F,[EE(k,4) EE(k,3)]),2)>=2);
        disp('ENTERED ELSE');
    end
    
    if (EE(k,2)>size(V_fine,1) && EE(k,4)>size(V_fine,1))
        [EE(k,2) EE(k,1)]
        [EE(k,4) EE(k,3)]
        disp('ENTERED CORRECTING NORMALS');
        if (EE(k,2)<EE(k,4))
            [EE(k,2) EE(k,1)]
            f = find(sum(ismember(F,[EE(k,2) EE(k,1)]),2)>=2);
        else
            [EE(k,4) EE(k,3)]
            f = find(sum(ismember(F,[EE(k,4) EE(k,3)]),2)>=2);
        end
%         f = find(sum(ismember(F,[EE(k,4) EE(k,3)]),2)>=2);
    end
    
    if (size(f,1)==2)
        N1 = normalizerow(normals([V0(F(f(1),1),:);V0(F(f(1),2),:);V0(F(f(1),3),:)],[1 2 3]));
        N2 = normalizerow(normals([V0(F(f(2),1),:);V0(F(f(2),2),:);V0(F(f(2),3),:)],[1 2 3]));
        N = normalizerow(N1+N2);
    else
        N1 = normals([V0(F(f(1),1),:);V0(F(f(1),2),:);V0(F(f(1),3),:)],[1 2 3]);
        N = N1/norm(N1);
    end
    
    xb = (1-u2)*A20+u2*B20;
    xa = (1-u1)*A10+u1*B10;

      disp('dot product at time 0');
      dot(xb-xa,N)
      
      xb = (1-u2)*A21+u2*B21;
      xa = (1-u1)*A11+u1*B11;

      disp('dot product at time 1');
      dot(xb-xa,N)
%       disp('distances at time 0 - edge-egde')
%       max(norm(P-xa),norm(P-xb))
      if (dot(xb-xa,N)<0 && ((EE(k,2)>size(V_fine,1) && EE(k,4)>size(V_fine,1))))
          disp('WARNING: dot(xb-xa,N)<0 and both belong to coarse mesh. Flipping normals')
          N = -N;
      end
%       if (dot(xb-xa,N)>0 && j<i && EE(k,4) <= size(V_fine,1))
%           disp('WARNING: dot(xb-xa,N)>0 and one of the vertices belong to fine mesh')
%           N = -N;
%       end
    
    EdgeColNX(i,j) = N(1);
    EdgeColNY(i,j) = N(2);
    EdgeColNZ(i,j) = N(3);

    
  end
  
%     % Initialize output collisions as zero
%   EdgeCol = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColX = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColY = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColZ = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColNX = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColNY = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColNZ = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColU1 = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColU2 = sparse(size(EE_all,1),size(EE_all,1));
%   EdgeColT = sparse(size(EE_all,1),size(EE_all,1));