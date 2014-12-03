function [V,F] = mcfskel_simplified(V0,F0,varargin)
  % MCFSKEL_SIMPLIFIED
  % mcfskel_simplified(V0,F0,varagin)
  %
  % Simpler version of "Mean Curvature Skeleletons"
  % by Tagliasacchi et al. (SGP 2012). Here we don't perform mesh 
  % operations and we also give the possibily of using cMCF
  %
  % Inputs:
  %   V0  list of surface vertex positions of exterior mesh, # vertices by 3
  %   F0  list of surface face indices of exterior triangle mesh, # faces by 3
  %     'delta_t' followed by timestep of the flow
  %     'pole_weight' weight for the contraction towards the poles
  %     'steps' number of steps of the flow.
  %     'lap_matrix' a matrix that will be used every iteration as the
  %     laplacian. If it is the cotmatrix of some initial mesh, then
  %     we are applying conformalized MCF instead of standard MCF.
  % Outputs:
  %   V  list of surface vertices
  %   F  list of faces (different from the input one)
  
  delta_t = 1e-4;
  pole_weight = 5e-3;
  steps = 100;
  lap_matrix = [];
  
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'delta_t'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              delta_t = varargin{ii};
          case 'pole_weight'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              pole_weight = varargin{ii};
          case 'steps'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              steps = varargin{ii};
          case 'lap_matrix'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              lap_matrix = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
  
  % decide which flow to use based on the input Laplacian matrix
  std_flow = true;
  if ~isempty(lap_matrix)
      std_flow = false;
      L = lap_matrix;
  end
  
  % initialize solution
  V = V0;
  F = F0;
  
  % calculate Voronoi poles
  [V_poles F_poles] = voronoi_poles_mcfskel(V0,F0);
  
  % loop 
  for k=1:steps
      
      k
      
      % re-compute Laplacian matrix
      if (std_flow)
          L = cotmatrix(V,F);
      end
      D = massmatrix(V,F,'full');
      
      % define system matrix (concatenation of MCF matrix and identity)
      S = [D-delta_t*L;pole_weight*speye(size(V0,1))];
      % augmented RHS
      V2 = [D*V;pole_weight*V_poles];
      
      % solve linear system
      V = S\(V2);
      
  end