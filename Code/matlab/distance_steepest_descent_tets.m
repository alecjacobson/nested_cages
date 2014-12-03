function [V,vertex_distances,moving_vertices] = distance_steepest_descent_tets(V0,F0,SD,DV,DT,dt,steps,w_lap,L,order,scheme,V_coarse,F_coarse,gradient,moving_vertices)

  % write description
  
  V_shrink = V0;
  n_fine = size(V_shrink,1);
    
  grad_vertices = zeros(size(V0));
  grad_energy = zeros(size(V0));
  vertex_distances = -ones(size(V0,1),1); 
  
  for k=1:steps
            
      % do intersection test to control the flow
      IF = intersect_other(V_shrink,F0,V_coarse,F_coarse);
      
      if (order==1)
          
          area_faces = doublearea(V_shrink,F0)/2;
          
          moving_vertices_idx = find(moving_vertices);
          
          tic 
          I = in_element_aabb(DV,DT,V_shrink(moving_vertices_idx,:));
          
          [II,~,IV] = find(I);
          V_cur = V_shrink(moving_vertices_idx,:);
          B = barycentric_coordinates(V_cur(II,:),DV(DT(IV,1),:),DV(DT(IV,2),:),DV(DT(IV,3),:),DV(DT(IV,4),:));
          B = sparse(repmat(II,1,size(B,2)),repmat(1:size(B,2),size(B,1),1),B,numel(I),size(B,2));
          toc

          grad_vertices(moving_vertices_idx,1) = gradient(I);
          grad_vertices(moving_vertices_idx,2) = gradient(size(DT,1)+I);
          grad_vertices(moving_vertices_idx,3) = gradient(2*size(DT,1)+I);
                    
          vertex_distances(moving_vertices_idx) = B(:,1).*SD(DT(I,1)) + B(:,2).*SD(DT(I,2)) + B(:,3).*SD(DT(I,3)) + B(:,4).*SD(DT(I,4));
                    
          disp('passou aqui 2.9');
          
          for i=moving_vertices_idx(1):moving_vertices_idx(end)
              
              face_idx = mod(find(F0==i)-1,size(F0,1))+1;
              % find adjacent vertices (includes the vertex itself)
              adj_vertices = unique(F0(face_idx,:));
              
              if (sum(vertex_distances(adj_vertices)<=-1.0)==size(adj_vertices,1) && sum(ismember(face_idx,IF(:,1)))==0)
                  
                  grad_energy(i,:) = 0;
                  moving_vertices(i) = 0;
                  
              else
                  
                  area_weight = sum(area_faces(face_idx))/size(face_idx,1);
                  grad_energy(i,:) = grad_vertices(i,:)*area_weight;
                  
              end
              
          end
                    
          disp('passou aqui 2.10');
          
      elseif (order==2)
          
          tic
          
          area_faces = doublearea(V_shrink,F0)/2;
          
          moving_vertices_idx = find(moving_vertices);
          moving_faces_idx = mod(find(ismember(F0,moving_vertices_idx))-1,size(F0,1))+1;
          
          % initialize functions for the face as zero
          dD_dx = zeros(size(F0,1),1);
          dD_dy = zeros(size(F0,1),1);
          dD_dz = zeros(size(F0,1),1);
          
          grad_vertex_x = zeros(size(F0,1),1);
          grad_vertex_y = zeros(size(F0,1),1);
          grad_vertex_z = zeros(size(F0,1),1);
          
          for j=1:3
          
              if (j==1)
                  % first quadrature point
                  p1 = 0.5*(V0(F0(moving_faces_idx,1),:)+V0(F0(moving_faces_idx,2),:));
              elseif (j==2)
                  % second quadrature point
                  p1 = 0.5*(V0(F0(moving_faces_idx,2),:)+V0(F0(moving_faces_idx,3),:));
              elseif (j==3)
                  % third quadrature point
                  p1 = 0.5*(V0(F0(moving_faces_idx,3),:)+V0(F0(moving_faces_idx,1),:));
              end
              
              I = in_element_aabb(DV,DT,p1);
              
              grad_vertex_x(moving_faces_idx) = gradient(I);
              grad_vertex_y(moving_faces_idx) = gradient(size(DT,1)+I);
              grad_vertex_z(moving_faces_idx) = gradient(2*size(DT,1)+I);

              dD_dx = dD_dx+(1/3)*grad_vertex_x;
              dD_dy = dD_dy+(1/3)*grad_vertex_y;
              dD_dz = dD_dz+(1/3)*grad_vertex_z;
              
          end
          
          grad_faces = [area_faces area_faces area_faces].*[dD_dx dD_dy dD_dz];
          
          % distance for vertices
          I = in_element_aabb(DV,DT,V_shrink(moving_vertices_idx,:));
          
          [II,~,IV] = find(I);
          V_cur = V_shrink(moving_vertices_idx,:);
          B = barycentric_coordinates(V_cur(II,:),DV(DT(IV,1),:),DV(DT(IV,2),:),DV(DT(IV,3),:),DV(DT(IV,4),:));
          B = sparse(repmat(II,1,size(B,2)),repmat(1:size(B,2),size(B,1),1),B,numel(I),size(B,2));
          
          vertex_distances(moving_vertices_idx) = B(:,1).*SD(DT(I,1)) + B(:,2).*SD(DT(I,2)) + B(:,3).*SD(DT(I,3)) + B(:,4).*SD(DT(I,4));
%           sum(max([SD(DT(I,1)) SD(DT(I,2)) SD(DT(I,3)) SD(DT(I,4))],[],2)==0)
%           vertex_distances(moving_vertices_idx) = max([SD(DT(I,1)) SD(DT(I,2)) SD(DT(I,3)) SD(DT(I,4))],[],2);
          
          for i=moving_vertices_idx(1):moving_vertices_idx(end)
              
              face_idx = mod(find(F0==i)-1,size(F0,1))+1;
              adj_vertices = unique(F0(face_idx,:));
                            
              if (sum(vertex_distances(adj_vertices)<-1.0)==size(adj_vertices,1) && sum(ismember(face_idx,IF(:,1)))==0)
    
                  grad_energy(i,:) = 0.0;
                  moving_vertices(i) = 0;
                  
              else
                  
                  grad_energy(i,:) = sum(grad_faces(face_idx,:).*(area_faces(face_idx)*ones(1,3)))/sum(area_faces(face_idx));
                  
              end
              
          end
          
          toc
          
      elseif (order==3)
          
          tic
          
          area_faces = doublearea(V_shrink,F0)/2;
          
          moving_vertices_idx = find(moving_vertices);
          moving_faces_idx = mod(find(ismember(F0,moving_vertices_idx))-1,size(F0,1))+1;
          
          % initialize functions for the face as zero
          dD_dx = zeros(size(F0,1),1);
          dD_dy = zeros(size(F0,1),1);
          dD_dz = zeros(size(F0,1),1);
          
          grad_vertex_x = zeros(size(F0,1),1);
          grad_vertex_y = zeros(size(F0,1),1);
          grad_vertex_z = zeros(size(F0,1),1);
          
          for j=1:4
                        
              if (j==1)
                  % first quadrature point
                  p1 = (1/3)*(V0(F0(moving_faces_idx,1),:)+V0(F0(moving_faces_idx,2),:)+V0(F0(moving_faces_idx,3),:));
              elseif (j==2)
                  % second quadrature point
                  p1 = (2/15)*V0(F0(moving_faces_idx,1),:)+(11/15)*V0(F0(moving_faces_idx,2),:)+(2/15)*V0(F0(moving_faces_idx,3),:);
              elseif (j==3)
                  % third quadrature point
                  p1 = (2/15)*V0(F0(moving_faces_idx,1),:)+(2/15)*V0(F0(moving_faces_idx,2),:)+(11/15)*V0(F0(moving_faces_idx,3),:);
              elseif (j==4)
                  % third quadrature point
                  p1 = (11/15)*V0(F0(moving_faces_idx,1),:)+(2/15)*V0(F0(moving_faces_idx,2),:)+(2/15)*V0(F0(moving_faces_idx,3),:);
              end

              I = in_element_aabb(DV,DT,p1);
              
              grad_vertex_x(moving_faces_idx) = gradient(I);
              grad_vertex_y(moving_faces_idx) = gradient(size(DT,1)+I);
              grad_vertex_z(moving_faces_idx) = gradient(2*size(DT,1)+I);
              
              if (j==1)
                  
                  dD_dx = dD_dx-(27/48)*grad_vertex_x;
                  dD_dy = dD_dy-(27/48)*grad_vertex_y;
                  dD_dz = dD_dz-(27/48)*grad_vertex_z;
                  
              else
                  
                  dD_dx = dD_dx+(25/48)*grad_vertex_x;
                  dD_dy = dD_dy+(25/48)*grad_vertex_y;
                  dD_dz = dD_dz+(25/48)*grad_vertex_z;
                  
              end
              
          end
          
          grad_faces = [area_faces area_faces area_faces].*[dD_dx dD_dy dD_dz];

          
          % distance for vertices
          I = in_element_aabb(DV,DT,V_shrink(moving_vertices_idx,:));
          [II,~,IV] = find(I);
          V_cur = V_shrink(moving_vertices_idx,:);
          B = barycentric_coordinates(V_cur(II,:),DV(DT(IV,1),:),DV(DT(IV,2),:),DV(DT(IV,3),:),DV(DT(IV,4),:));
          B = sparse(repmat(II,1,size(B,2)),repmat(1:size(B,2),size(B,1),1),B,numel(I),size(B,2));
          vertex_distances(moving_vertices_idx) = B(:,1).*SD(DT(I,1)) + B(:,2).*SD(DT(I,2)) + B(:,3).*SD(DT(I,3)) + B(:,4).*SD(DT(I,4));
          
          for i=moving_vertices_idx(1):moving_vertices_idx(end)
              
              face_idx = mod(find(F0==i)-1,size(F0,1))+1;
              adj_vertices = unique(F0(face_idx,:));
              
              if (sum(vertex_distances(adj_vertices)<-1.0)==size(adj_vertices,1) && sum(ismember(face_idx,IF(:,1)))==0)
                  
                  grad_energy(i,:) = 0.0;
                  moving_vertices(i) = 0;
                  
              else
                  
                  grad_energy(i,:) = sum(grad_faces(face_idx,:).*(area_faces(face_idx)*ones(1,3)))/sum(area_faces(face_idx));
                  
              end
              
          end
          
          toc
          
          
      else
          
          error('integration order between 1 and 3 are available')
          
              
      end
      
      M = massmatrix(V_shrink,F0,'voronoi');
      if (strcmp(scheme,'explicit'))
          
          grad_lap = -M\(L*V_shrink);
          V_shrink = V_shrink - dt*(grad_energy+w_lap*grad_lap);
          
      elseif (strcmp(scheme,'semi-implicit'))
          
%           moving_vertices_idx = find(moving_vertices);
          
          S = M-dt*w_lap*L;
          V_shrink_solve = S\(M*V_shrink -dt*(grad_energy));
          
%           V_shrink(moving_vertices_idx,:) = V_shrink_solve(moving_vertices_idx,:);
          V_shrink = V_shrink_solve;
          
      else
          
          error('unknown integration scheme')
          
      end
      
      disp('passou aqui 2.11');
      
      
  end
  
  V = V_shrink;
      
      
      
      
      
      
%   % old version (without quadrature) 
%   % distance for points on the surface (bilinear interpolation)
%   [I,B1,B2,B3] = in_element(DV,DT,V_shirnk);  
%   V_SD = B1*SD(DT(:,1),:) + B2*SD(DT(:,2),:) + B3*SD(DT(:,3),:) +...
%       (I-B1-B2-B3)*SD(DT(:,4),:);
%   % distance for points on the surface displaced on the x-dir by +delta_der 
%   [I,B1,B2,B3] = in_element(DV,DT,V_shirnk+dx*[ones(n_fine,1) zeros(n_fine,1) zeros(n_fine,1)]);  
%   V_SD_dx = B1*SD(DT(:,1),:) + B2*SD(DT(:,2),:) + B3*SD(DT(:,3),:) +...
%       (I-B1-B2-B3)*SD(DT(:,4),:);
%   % distance for points on the surface displaced on the x-dir by -delta_der 
%   [I,B1,B2,B3] = in_element(DV,DT,V_shirnk-dx*[ones(n_fine,1) zeros(n_fine,1) zeros(n_fine,1)]);  
%   V_SD_dx_minus = B1*SD(DT(:,1),:) + B2*SD(DT(:,2),:) + B3*SD(DT(:,3),:) +...
%       (I-B1-B2-B3)*SD(DT(:,4),:);