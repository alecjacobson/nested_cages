function [V,moving_vertices,grad_V,quad_points,grad_quadrature,C,s] = signed_distance_direction_quadrature_matrix(V0,F0,steps,order,V_coarse,F_coarse,moving_vertices,A,M,w_lap,L,sign_increase,varargin)
  % SIGNED_DISTANCE_DIRECTION_QUADRATURE_MATRIX
  % [V,moving_vertices,grad_V,quad_points,grad_quadrature,C] = signed_distance_direction_quadrature_matrix(V0,F0,steps,order,V_coarse,F_coarse,moving_vertices,A,M,w_lap,L)
  %
  % Given an initial triangle mesh (V0,F0) and a simplified mesh 
  % (V_coarse,F_coarse), it shrinks (V0,F0) using the signed distance
  % direction as the gradient plus line search. This version is faster
  % due to avoiding for loops and doing everything on matrix form
  %
  % Input:
  %   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
  %   F0  (#faces)x3 list of vertex indices that form each face of the
  %   initial mesh
  %   steps: number of steps to take
  %   order: quadrature order for the intergaryion (order = 1, 2 or 3)
  %   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
  %   coarse mesh
  %   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
  %   of the coarse mesh
  %   moving_vertices: (#vertices)x1 vector of bools specifing which
  %   vertices should move
  %   'A': matrix that converts gradients at quadrature points to gradients
  %   at fine mesh vertices
  %   M: massmatrix used for stepping the flow (usually defined w.r.t.
  %   initial mesh)
  %   w_lap: scalar that defines the ammount of smoothing
  %   L: cotangent matrix w.r.t. initial mesh
  %   sign_increase: 1 of maximum decrase and -1 if maximum increase
  % Output:
  %   V   (#vertices)x3 list of new mesh vertex positions of 
  %   the fine mesh
  %   moving_vertices: (#vertices)x1 vector of bools specifing which
  %   vertices should move (updated)
  %   grad_V: (#vertices)x3 list of gradients at mesh vertices
  %   quad_points: list of quadrature points used for integration
  %          if quadrature_order == 1, then size(quad_points) = [#faces,3]
  %          if quadrature_order == 2, then size(quad_points) = [3*#faces,3]
  %          if quadrature_order == 3, then size(quad_points) = [4*#faces,3]
  %   grad_quadrature: list of gradients at quadrature points
  %          if quadrature_order == 1, then size(grad_quadrature) = [#faces,3]
  %          if quadrature_order == 2, then size(grad_quadrature) = [3*#faces,3]
  %          if quadrature_order == 3, then size(grad_quadrature) = [4*#faces,3]
  %   C: closest points on the coarse mesh to quadrature points (same size
  %      as 'quad_points'
  %   s: resulting step size
  
  V_shrink = V0;
  
  % initialize scalar factor for line search
  s = 1e-3;
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'step'
              ii = ii+1;
              assert(ii<=numel(varargin));
              s = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  initial_s = s;
  
  % calculate traingles that are not totally inside
  moving_faces = (1:size(F0,1))';
  
  for k=1:steps
  
      if (order==1)
          
          % only quadrature point (barycenters)
          p123 = (1/3)*(V0(F0(moving_faces,1),:)+V0(F0(moving_faces,2),:)+V0(F0(moving_faces,3),:));
          quad_points = p123;
          
          [grad_p_all,~,C] = signed_distance_direction(quad_points,V_coarse,F_coarse);
          grad_quadrature = -grad_p_all;
             
          A_not_in = A(:,moving_faces);
          grad_V = (M\(A_not_in*grad_quadrature));
         
      elseif (order==2)
          
          % define the 3 quadrature points (edge-midpoints)
          p12 = 0.5*(V0(F0(moving_faces,1),:)+V0(F0(moving_faces,2),:));
          p23 = 0.5*(V0(F0(moving_faces,2),:)+V0(F0(moving_faces,3),:));
          p31 = 0.5*(V0(F0(moving_faces,3),:)+V0(F0(moving_faces,1),:));
          quad_points = [p12;p23;p31];
          
          [grad_p_all,~,C] = signed_distance_direction(quad_points,V_coarse,F_coarse);
          grad_quadrature = -grad_p_all;
          
          A_not_in = A(:,[moving_faces;size(F0,1)+moving_faces;2*size(F0,1)+moving_faces]);
          grad_V = (M\(A_not_in*grad_quadrature));
          
      elseif (order==3)
          
          % define the 4 quadrature points
          p1 = (1/3)*(V0(F0(moving_faces,1),:)+V0(F0(moving_faces,2),:)+V0(F0(moving_faces,3),:));
          p2 = (2/15)*V0(F0(moving_faces,1),:)+(11/15)*V0(F0(moving_faces,2),:)+(2/15)*V0(F0(moving_faces,3),:);
          p3 = (2/15)*V0(F0(moving_faces,1),:)+(2/15)*V0(F0(moving_faces,2),:)+(11/15)*V0(F0(moving_faces,3),:);
          p4 = (11/15)*V0(F0(moving_faces,1),:)+(2/15)*V0(F0(moving_faces,2),:)+(2/15)*V0(F0(moving_faces,3),:);
          quad_points = [p1;p2;p3;p4];
          
          [grad_p_all,~,C] = signed_distance_direction(quad_points,V_coarse,F_coarse);
          grad_quadrature = -grad_p_all;
          
          A_not_in = A(:,[moving_faces;size(F0,1)+moving_faces;2*size(F0,1)+moving_faces;3*size(F0,1)+moving_faces]);
          grad_V = (M\(A_not_in*grad_quadrature));
          
      else
          
          error('integration order between 1 and 3 are available')
          
      end
            
      % smooth if necessary (changed by constant smoothing)
      S = s*ones(size(V0,1),1)*sign_increase;
      M_new = massmatrix(V0,F0,'barycentric');
      rhs_MCF = M_new*V_shrink;
      A_MCF = (M_new-initial_s*L);
      V_shrink = [M*speye(size(V0,1));w_lap*A_MCF]\[M*(V_shrink-(S*ones(1,3)).*(grad_V));w_lap*rhs_MCF];
      
  end
  
  % output gradients times steps
  grad_V = (S*ones(1,3)).*(grad_V);
  V = V_shrink;
