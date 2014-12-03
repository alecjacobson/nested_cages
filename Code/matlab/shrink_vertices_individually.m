function [V_shrink] = shrink_vertices_individually(V0,F0,V_coarse,F_coarse,moving_vertices)
  % SHRINK_VERTICES_INDIVIDUALLY
  % V_shrink = shrink_vertices_individually(V0,F0,V_coarse,F_coarse,moving_vertices)
  
  % Given some vertices of the mesh, shrink them individually moving
  % them towards their gradient directions
  
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
  %   coarse mesh
  %   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
  %   of the coarse mesh
  %   moving_vertices: (#moving_vertices)x1 vector with indices of the
  %   vertices that should move
  % Input:
  %   V_shrink  (#vertices)x3 list of mesh vertex positions of the final fine mesh
  
  V_shrink = V0;
  
  % initialize scalar factor for line search
  initial_s = 0.005;
  s = initial_s;
  % factor to impose for local decrease (Armijo)
  c = 1e-7;
  factor = 0.8;
  
  N = normalsvertex(V0,F0);
  Grad = zeros(size(moving_vertices,2),3);
  
  for k=1:size(moving_vertices,2)
      
      p = V_shrink(moving_vertices(k),:);
      
      
          [Q,~] = qr(N(moving_vertices(k),:)');
          % centered mesh (old coordinates)
          centered_mesh = V_coarse-ones(size(V_coarse,1),1)*p;
          % coordinates on the basis Q
          centered_mesh = (Q*diag([0.5 1 1])*Q'*centered_mesh')';
          
          [SD,~,C] = signed_distance(p-p,centered_mesh,F_coarse);
          C = (Q'*diag([2 1 1])*Q*C')';
  

%           [SD,~,C] = signed_distance(p,V_coarse,F_coarse);

          if (SD>0)
              V_shrink(moving_vertices(k),:) = p+C;
          end
          
      
%       grad = -signed_distance_direction(p,V_coarse,F_coarse);
      
      
%       local_slope = -sum(sum(grad.*grad));
%       
%       while (true)
%           dir_attempt = -s*(grad);
%           p_attempt = p + dir_attempt;
%           
%           [SD_attempt,~,~] = signed_distance(p_attempt,V_coarse,F_coarse);
%           
%           if (SD_attempt>SD+c*s*local_slope)
%               s = factor*s;
%           else
%               s = initial_s;
%               break;
%           end
%           
%       end
%       
%       V_shrink(moving_vertices(k),:) = V0(moving_vertices(k),:)-s*grad;
      
  end
  
