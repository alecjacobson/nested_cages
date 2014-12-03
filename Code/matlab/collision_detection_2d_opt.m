function [Col,ColX,ColY,ColNX,ColNY,ColU,ColT] = collision_detection_2d_opt(V0,V1,E,varargin)
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
  %   Optional:
  %     'last_only' test only if last point of the mesh collides
  % Output:
  %   Col   MxN matrix such that Col(i,j) == 1 iff there is a collision between
  %     vertex i and edges j
  %   ColX   MxN matrix such that ColY(i,j) is the y-coordinates
  %     2D position where vertex i collides with edge j. 
  %   ColY   MxN matrix such that ColY(i,j) is the y-coordinates
  %     2D position where vertex i collides with edge j. 
  %   ColNX   MxN matrix such that ColNX(i,j) is the x-coordinate
  %     of the collision normal of vertex i collides with edge j. 
  %   ColNY   MxN matrix such that ColNY(i,j) is the y-coordinate
  %     of the collision normal of vertex i collides with edge j. 
  %   ColU   MxN matrix such that ColU(i,j) is the fraction of edge
  %      j where vertex i collides. 
  %   ColT   MxN matrix such that ColT(i,j) is the time instant when
  %   vertex i collides with edge j.
  
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
  ColNX = sparse(size(V0,1),size(E,2));
  ColNY = sparse(size(V0,1),size(E,2));
  ColU = sparse(size(V0,1),size(E,2));
  ColT = sparse(size(V0,1),size(E,2));
  % Initialize coefficients of the polynomial to be solved as zeros
  p = zeros(1,3);

  % number of edges
  m=size(E,2);
  % number of vertices
  n = size(V0,1);
  % indices of start of groups
  bid = round(linspace(1,m+1,(m/2)+1));
  %bid = [1 m+1];
  % number of bounding boxes
  nb = numel(bid)-1;
  % Bounding boxes: B(b,:) = [xmin xmax ymin ymax] for box b
  B = zeros(nb,4);
  % B2V(b,:) reveals vertices in bounding box v
  B2V = sparse(nb,n);
  for b = 1:nb
    % assumes all vertices have at least one incident edge!
    v = reshape(E(:,bid(b):bid(b+1)-1),2*(bid(b+1)-bid(b)),1);
    B(b,1) = min([V0(v,1);V1(v,1)]);
    B(b,2) = max([V0(v,1);V1(v,1)]);
    B(b,3) = min([V0(v,2);V1(v,2)]);
    B(b,4) = max([V0(v,2);V1(v,2)]);
    B2V(b,v) = 1;
  end


  % loop over this bounding boxes
  for b1 = 1:nb
      
    % loop over that bounding boxes
    for b2 = b1:nb
        
      if b1 ~= b2
        % if this's min is greater than that's max or that's min is greater
        % than this's max then we can skip it
        if ...
            B(b1,1) > B(b2,2) || ...
            B(b2,1) > B(b1,2) || ...
            B(b1,3) > B(b2,4) || ...
            B(b2,3) > B(b1,4)
          continue;
        end 
        % edges in either box
        J = [bid(b1):bid(b1+1)-1,bid(b2):bid(b2+1)-1];
      else
        J = [bid(b1):bid(b1+1)-1];
      end
      
      % Recover original edge index
      EE = E(:,J);
      % vertices in either box
      IV = find(B2V(b1,:)|B2V(b2,:));
  
      % loop over all edges j and vertices i
      for j=1:size(EE,2)
          E1j = EE(1,j);
          E2j = EE(2,j);
          
          for i=IV
            % faster than setdiff
            if i == E1j || i == E2j
              continue;
            end
            
            % go to next iteration if vertex is not the last one, and 
            % it is set to test only the last one
            if last_only==1 && i~=size(V0,1)
                continue;
            end
              
              % precompute inner products that will be needed
              inner00 = sum(Normals_0(J(j),:).*(V0(i,:)-V0(E1j,:)));
              inner01 = sum(Normals_0(J(j),:).*(V1(i,:)-V1(E1j,:)));
              inner10 = sum(Normals_1(J(j),:).*(V0(i,:)-V0(E1j,:)));
              inner11 = sum(Normals_1(J(j),:).*(V1(i,:)-V1(E1j,:)));
              
              % calculate coefficients of the polynomial
              p(1) = inner00-inner01-inner10+inner11;
              p(2) = -2*inner00+inner01+inner10;
              p(3) = inner00;
              
              if p(1) == 0 && p(2) == 0
                continue;
              elseif p(1) == 0
                sols = -p(3)/p(2);
              else
                sols = zeros(2,1);
                % http://en.wikipedia.org/wiki/Quadratic_equation#Floating-point_implementation
                d = sqrt(p(2)^2 - 4*p(1)*p(3));
                sols(1) = (-p(2) - sign(p(2))*d)/(2*p(1));
                sols(2) = p(3)/(p(1)*sols(1));
              end
              
              % select roots that belong to [0,1]. If more than one,
              % select the smallest one. Also exclude non real solutions
              t = sols((sols>0.0)&(sols<1.0)&imag(sols)==0);
              t = min(t);
              
              % if returned some t from previous selection
              if(size(t,1)==1)
                                    
                  % compute all necessary points
                  P_t = (1-t)*V0(i,:)+t*V1(i,:);
                  A_t = (1-t)*V0(E1j,:)+t*V1(E1j,:);
                  B_t = (1-t)*V0(E2j,:)+t*V1(E2j,:);

                  % determine if P lies on AB edge (it does if u below belongs
                  % to [0,1])
                  [~,idx] = max(abs(B_t-A_t));
                  u = (P_t(idx)-A_t(idx))/(B_t(idx)-A_t(idx));

                  % if it u belongs to [0,1], output colliding point,
                  % normal and fraction of the edge
                  
                  if (u>=0&&u<=1)
                      
                       Col(i, J(j)) = 1;
                       ColX(i,J(j)) = P_t(1);
                       ColY(i,J(j)) = P_t(2);
                       
%                        ColNX(i,J(j)) = Normals_1(J(j),1);
%                        ColNY(i,J(j)) = Normals_1(J(j),2);
%                        norm_vec = norm(Normals_1(J(j),:));


                       if (i==1)
                           ColNX(i,J(j)) = Normals_1(n,1)+Normals_1(i,1);
                           ColNY(i,J(j)) = Normals_1(n,2)+Normals_1(i,2);
                           norm_vec = norm(Normals_1(n,:)+Normals_1(i,:));
                       else
                           ColNX(i,J(j)) = Normals_1(i-1,1)+Normals_1(i,1);
                           ColNY(i,J(j)) = Normals_1(i-1,2)+Normals_1(i,2);
                           norm_vec = norm(Normals_1(i-1,:)+Normals_1(i,:));
                       end
                       
                       ColNX(i,J(j)) = ColNX(i,J(j))/norm_vec;
                       ColNY(i,J(j)) = ColNY(i,J(j))/norm_vec;
                       
%                        % CCW test on the initial points to
%                        % determine if collision normal has to be flipped
%                        ccw_test = det([V0(i,1) V0(E1j,1) V0(E2j,1);...
%                            V0(i,2) V0(E1j,2) V0(E2j,2);1 1 1])>=0;
%                        if (ccw_test)
% %                             ColNX(i,J(j)) = (1-t)*Normals_0(J(j),1)+t*Normals_1(J(j),1);
% %                             ColNY(i,J(j)) = (1-t)*Normals_0(J(j),2)+t*Normals_1(J(j),2);
% %                             norm_vec = norm((1-t)*Normals_0(J(j),:)+t*Normals_1(J(j),:));
% %                             ColNX(i,J(j)) = ColNX(i,J(j))/norm_vec;
% %                             ColNY(i,J(j)) = ColNY(i,J(j))/norm_vec;
% 
%                             % normals have to be the ones at the end of the
%                             % time step
%                             ColNX(i,J(j)) = Normals_1(J(j),1);
%                             ColNY(i,J(j)) = Normals_1(J(j),2);
%                             norm_vec = norm(Normals_1(J(j),:));
%                             ColNX(i,J(j)) = ColNX(i,J(j))/norm_vec;
%                             ColNY(i,J(j)) = ColNY(i,J(j))/norm_vec;
%                        else
% %                             ColNX(i,J(j)) = -((1-t)*Normals_0(J(j),1)+t*Normals_1(J(j),1));
% %                             ColNY(i,J(j)) = -((1-t)*Normals_0(J(j),2)+t*Normals_1(J(j),2));
% %                             norm_vec = norm(-((1-t)*Normals_0(J(j),:)+t*Normals_1(J(j),:)));
% %                             ColNX(i,J(j)) = ColNX(i,J(j))/norm_vec;
% %                             ColNY(i,J(j)) = ColNY(i,J(j))/norm_vec;
% 
%                             disp('entered here!');
%                             % normals have to be the ones at the end of the
%                             % time step
%                             ColNX(i,J(j)) = -Normals_1(J(j),1);
%                             ColNY(i,J(j)) = -Normals_1(J(j),2);
%                             norm_vec = norm(-Normals_1(J(j),:));
%                             ColNX(i,J(j)) = ColNX(i,J(j))/norm_vec;
%                             ColNY(i,J(j)) = ColNY(i,J(j))/norm_vec;
%                        end

                       ColU(i,J(j)) = u;
                       ColT(i,J(j)) = t;
                  end
                  
              end
              
          end
          
        end
    end
  end
