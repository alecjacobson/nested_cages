function A = grad_quadrature_to_vertices(V0,F0,area_initial,quadrature_order)
  % GRAD_QUADRATURE_TO_VERTICES
  % grad_quadrature_to_vertices(V0,F0,areas,quadrature_order)
  %
  % Given a triangle mesh (V0,F0), associated triangle 'areas' and
  % a 'quadrature_order', defines the mesh that transform gradient at
  % face quadrature vertives to mesh vertices.
  % 
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   area_initial: (#faces)x1 triangle ares
  %   quadrature_order: quadrature order for the intergaryion 
  %   (quadrature_order = 1, 2 or 3)
  % Output:
  %   A: if (quadrature_order==1) 
  %        (#vertices)x(#faces) that converts gradient at barycenters to 
  %        vertices on the mesh
  %      elseif (quadrature_order==2)
  %        (#vertices)x(3*#faces) that converts gradient at 2nd order 
  %        quadrature points to vertices on the mesh
  %      elseif (quadrature_order==3)
  %        (#vertices)x(4*#faces) that converts gradient at 3rd order 
  %        quadrature points to vertices on the mesh
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Check for cached result, do NOT edit variables until cache is checked,
  % your function code comes later. See below
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [cache_exists,cache_name] = find_cache();
  if cache_exists
    fprintf('Using cache...\n');
    load(cache_name);
    return;
  end
  fprintf('First time. Creating cache...\n');
  
  if (quadrature_order==1)
      % initialize matrix A (this should be in a separate function)
%       A = sparse(size(V0,1),size(F0,1));
      
      I = zeros(10*size(V0,1),1);
      J = zeros(10*size(V0,1),1);
      Val = zeros(10*size(V0,1),1);
      nnz_cur = 0;
      
      % defines A entries
      for i=1:size(V0,1)
          
          face_idx = mod(find(F0==i)-1,size(F0,1))+1;
%           A(i,face_idx) = area_initial(face_idx)/3;
          
          num_entries = size(face_idx,1);
          I(nnz_cur+1:nnz_cur+num_entries) = i;
          J(nnz_cur+1:nnz_cur+num_entries) = face_idx;
          Val(nnz_cur+1:nnz_cur+num_entries) = area_initial(face_idx)/3;
          
          nnz_cur = nnz_cur + num_entries;
          
      end
      
      I = nonzeros(I);
      J = nonzeros(J);
      Val = nonzeros(Val);
      A = sparse(I,J,Val);
            
  elseif (quadrature_order==2)
      
      % initialize matrix A (this should be in a separate function)
%       A = sparse(size(V0,1),3*size(F0,1));
      
      I = zeros(3*10*size(V0,1),1);
      J = zeros(3*10*size(V0,1),1);
      Val = zeros(3*10*size(V0,1),1);
      nnz_cur = 0;
      
      % defines A entries
      for i=1:size(V0,1)
          
          face_idx = mod(find(F0==i)-1,size(F0,1))+1;
          % find index of each vertex in the faces it belongs (to
          % select the write quadrature points)
          vtx_idx = mod(find((F0(face_idx,:))'==i)-1,3)+1;
          
          for m=1:size(vtx_idx,1)
              
              if (vtx_idx(m)==1)
%                   A(i,face_idx(m)) = A(i,face_idx(m))+area_initial(face_idx(m))/6;
%                   A(i,2*size(F0,1)+face_idx(m)) = A(i,2*size(F0,1)+face_idx(m))+area_initial(face_idx(m))/6;
                  
                  I(nnz_cur+1) = i;
                  J(nnz_cur+1) = face_idx(m);
                  Val(nnz_cur+1) = area_initial(face_idx(m))/6;
                  
                  I(nnz_cur+2) = i;
                  J(nnz_cur+2) = 2*size(F0,1)+face_idx(m);
                  Val(nnz_cur+2) = area_initial(face_idx(m))/6;
                  
              elseif (vtx_idx(m)==2)
%                   A(i,face_idx(m)) = A(i,face_idx(m)) + area_initial(face_idx(m))/6;
%                   A(i,size(F0,1)+face_idx(m)) = A(i,size(F0,1)+face_idx(m)) + area_initial(face_idx(m))/6;
                  
                  I(nnz_cur+1) = i;
                  J(nnz_cur+1) = face_idx(m);
                  Val(nnz_cur+1) = area_initial(face_idx(m))/6;
                  
                  I(nnz_cur+2) = i;
                  J(nnz_cur+2) = size(F0,1)+face_idx(m);
                  Val(nnz_cur+2) = area_initial(face_idx(m))/6;
                  
              elseif (vtx_idx(m)==3)
%                   A(i,size(F0,1)+face_idx(m)) = A(i,size(F0,1)+face_idx(m)) + area_initial(face_idx(m))/6;
%                   A(i,2*size(F0,1)+face_idx(m)) = A(i,2*size(F0,1)+face_idx(m)) + area_initial(face_idx(m))/6;
                  
                  I(nnz_cur+1) = i;
                  J(nnz_cur+1) = size(F0,1)+face_idx(m);
                  Val(nnz_cur+1) = area_initial(face_idx(m))/6;
                  
                  I(nnz_cur+2) = i;
                  J(nnz_cur+2) = 2*size(F0,1)+face_idx(m);
                  Val(nnz_cur+2) = area_initial(face_idx(m))/6;
                  
              else
                  
                  error('vertex index is not 1 nor 2 nor 3')
                  
              end
              
              nnz_cur = nnz_cur+2;
              
          end
          
      end
      
      I = nonzeros(I);
      J = nonzeros(J);
      Val = nonzeros(Val);
      A = sparse(I,J,Val);
      
      
  elseif (quadrature_order==3)
      
      % initialize matrix A
%       A = sparse(size(V0,1),4*size(F0,1));
      
      I = zeros(4*10*size(V0,1),1);
      J = zeros(4*10*size(V0,1),1);
      Val = zeros(4*10*size(V0,1),1);
      nnz_cur = 0;
      
      % defines A entries
      for i=1:size(V0,1)
          
          face_idx = mod(find(F0==i)-1,size(F0,1))+1;
          % find index of each vertex in the faces it belongs (to
          % select the write quadrature points)
          vtx_idx = mod(find((F0(face_idx,:))'==i)-1,3)+1;
          
          for m=1:size(vtx_idx,1)
              
              if (vtx_idx(m)==1)
%                   A(i,face_idx(m)) = A(i,face_idx(m)) + (-3/16)*area_initial(face_idx(m));
%                   A(i,size(F0,1)+face_idx(m)) = A(i,size(F0,1)+face_idx(m)) + (5/72)*area_initial(face_idx(m));
%                   A(i,2*size(F0,1)+face_idx(m)) = A(i,2*size(F0,1)+face_idx(m)) + (5/72)*area_initial(face_idx(m));
%                   A(i,3*size(F0,1)+face_idx(m)) = A(i,3*size(F0,1)+face_idx(m)) + (55/144)*area_initial(face_idx(m));
                  
                  I(nnz_cur+1) = i;
                  J(nnz_cur+1) = face_idx(m);
                  Val(nnz_cur+1) = (-3/16)*area_initial(face_idx(m));
                  
                  I(nnz_cur+2) = i;
                  J(nnz_cur+2) = size(F0,1)+face_idx(m);
                  Val(nnz_cur+2) = (5/72)*area_initial(face_idx(m));
                  
                  I(nnz_cur+3) = i;
                  J(nnz_cur+3) = 2*size(F0,1)+face_idx(m);
                  Val(nnz_cur+3) = (5/72)*area_initial(face_idx(m));
                  
                  I(nnz_cur+4) = i;
                  J(nnz_cur+4) = 3*size(F0,1)+face_idx(m);
                  Val(nnz_cur+4) = (55/144)*area_initial(face_idx(m));
                  
              elseif (vtx_idx(m)==2)
%                   A(i,face_idx(m)) = A(i,face_idx(m)) + (-3/16)*area_initial(face_idx(m));
%                   A(i,size(F0,1)+face_idx(m)) = A(i,size(F0,1)+face_idx(m)) + (55/144)*area_initial(face_idx(m));
%                   A(i,2*size(F0,1)+face_idx(m)) = A(i,2*size(F0,1)+face_idx(m)) + (5/72)*area_initial(face_idx(m));
%                   A(i,3*size(F0,1)+face_idx(m)) = A(i,3*size(F0,1)+face_idx(m)) + (5/72)*area_initial(face_idx(m));
                  
                  I(nnz_cur+1) = i;
                  J(nnz_cur+1) = face_idx(m);
                  Val(nnz_cur+1) = (-3/16)*area_initial(face_idx(m));
                  
                  I(nnz_cur+2) = i;
                  J(nnz_cur+2) = size(F0,1)+face_idx(m);
                  Val(nnz_cur+2) = (55/144)*area_initial(face_idx(m));
                  
                  I(nnz_cur+3) = i;
                  J(nnz_cur+3) = 2*size(F0,1)+face_idx(m);
                  Val(nnz_cur+3) = (5/72)*area_initial(face_idx(m));
                  
                  I(nnz_cur+4) = i;
                  J(nnz_cur+4) = 3*size(F0,1)+face_idx(m);
                  Val(nnz_cur+4) = (5/72)*area_initial(face_idx(m));
                  
              elseif (vtx_idx(m)==3)
%                   A(i,face_idx(m)) = A(i,face_idx(m)) + (-3/16)*area_initial(face_idx(m));
%                   A(i,size(F0,1)+face_idx(m)) = A(i,size(F0,1)+face_idx(m)) + (5/72)*area_initial(face_idx(m));
%                   A(i,2*size(F0,1)+face_idx(m)) = A(i,2*size(F0,1)+face_idx(m)) + (55/144)*area_initial(face_idx(m));
%                   A(i,3*size(F0,1)+face_idx(m)) = A(i,3*size(F0,1)+face_idx(m)) + (5/72)*area_initial(face_idx(m));
                  
                  I(nnz_cur+1) = i;
                  J(nnz_cur+1) = face_idx(m);
                  Val(nnz_cur+1) = (-3/16)*area_initial(face_idx(m));
                  
                  I(nnz_cur+2) = i;
                  J(nnz_cur+2) = size(F0,1)+face_idx(m);
                  Val(nnz_cur+2) = (5/72)*area_initial(face_idx(m));
                  
                  I(nnz_cur+3) = i;
                  J(nnz_cur+3) = 2*size(F0,1)+face_idx(m);
                  Val(nnz_cur+3) = (55/144)*area_initial(face_idx(m));
                  
                  I(nnz_cur+4) = i;
                  J(nnz_cur+4) = 3*size(F0,1)+face_idx(m);
                  Val(nnz_cur+4) = (5/72)*area_initial(face_idx(m));
                  
              else
                  
                  error('vertex index is not 1 nor 2 nor 3 nor 4')
                  
              end
              
              nnz_cur = nnz_cur+4;
              
          end
          
      end
      
      I = nonzeros(I);
      J = nonzeros(J);
      Val = nonzeros(Val);
      A = sparse(I,J,Val);
            
  else
      
      error('specify quadrature order 1, 2 or 3 ')
      
  end
  
  create_cache(cache_name);
  
end
  