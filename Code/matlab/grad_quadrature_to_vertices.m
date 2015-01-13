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

  if (quadrature_order==1)

      n = size(V0,1);
      m = size(F0,1);
      I = (1:m)';
      % see also signed_distance_direction_quadrature_matrix.m
      % p123 depends on F0(:,1), F0(:,2) and F0(:,3)
      
      A = ...
      sparse( ...
        [F0(:,1) F0(:,2) F0(:,3)], ...
        [I(:,1) I(:,1) I(:,1)], ...
        repmat(area_initial/3,1,3), ...
        n,m);

  elseif (quadrature_order==2)
    n = size(V0,1);
    m = size(F0,1);
    I = reshape((1:3*m)',m,3);
    % see also signed_distance_direction_quadrature_matrix.m
    % p12 and p31 depend on F0(:,1)
    % p23 and p12 depend on F0(:,2)
    % p31 and p23 depend on F0(:,3)
    A = ...
      sparse( ...
        [ F0(:,[1 1]) F0(:,[2 2]) F0(:,[3 3])], ...
        [  I(:,[1 3])  I(:,[2 1])  I(:,[3 2])], ...
        repmat(area_initial/6,1,6), ...
        n,3*m);

  elseif (quadrature_order==3)
      
    n = size(V0,1);
    m = size(F0,1);
    I = reshape((1:4*m)',m,4);
    % see also signed_distance_direction_quadrature_matrix.m
    % p1, p2, p3, p4 depend on F0(:,1)
    % p1, p2, p3, p4 depend on F0(:,2)
    % p1, p2, p3, p4 depend on F0(:,3)
    A = ...
      sparse( ...
        [ F0(:,[1 1 1 1]) F0(:,[2 2 2 2]) F0(:,[3 3 3 3])], ...
        [  I(:,[1 2 3 4])  I(:,[1 2 3 4])  I(:,[1 2 3 4])], ...
        repmat([-3/16 5/72 5/72 55/144 -3/16 55/144 5/72 5/72 -3/16 5/72 55/144 5/72],size(F0,1),1).*repmat(area_initial,1,12), ...
        n,4*m);
  else

      error('specify quadrature order 1, 2 or 3 ')

  end

end

