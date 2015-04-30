function Pall = cmcf_shrink(V,F,varargin)
  % CMCF_SHRINK Flow a surface according to "Can mean
  % curvature flow be made non-singular?" [Kazhdan et al. 2012]
  % Obs.: This function is a simple modification of 
  % gptoolbox/mesh/conformalized_mean_curvature_flow,
  % where we remove re-normalization (to actually shrink the mesh)
  % and output all states of the flow
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'MaxIter' followed by maximum number of iterations {100}
  %     'delta' followed by delta value, should roughly be in range
  %       [1e-13,1e13] {1}
  %     'LaplacianType' followed by 'cotangent' of 'uniform'.
  %     'V0' followed by #V by 3 mesh positions to treat as initial mesh (to
  %       build laplacian from)
  % Outputs:
  %   Pall  #V by 3 by #steps list of vertex positions
  %

  % default values
  delta = 1;
  max_iter = 100;
  laplacian_type = 'cotangent';
  V0 = V;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','delta','LaplacianType','V0'}, ...
    {'max_iter','delta','laplacian_type','V0'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  switch laplacian_type
  case 'cotangent'
    L = cotmatrix(V,F);
  case 'uniform'
    A = adjacency_matrix(F);
    L = A - diag(sparse(sum(A,2)));
  end

  U = V;
  Pall = zeros(size(V,1),3,max_iter+1);
  Pall(:,:,1) = U;
  iter = 1;
  while true
    U_prev = U;
    % 'full' seems slight more stable than 'barycentric' which is more stable
    % than 'voronoi'
    M = massmatrix(U,F,'barycentric');
    U = (M-delta*L)\(M*U);
    Pall(:,:,iter+1) = U;
    if iter >= max_iter
      break;
    end
    iter = iter + 1;
  end

end