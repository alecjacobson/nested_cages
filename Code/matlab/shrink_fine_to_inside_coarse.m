function Sall = shrink_fine_to_inside_coarse(P,P_coarse,varargin)
  % SHRINK_FINE_TO_INSIDE_COARSE
  % shrink_fine_to_inside_coarse(P,P_coarse,varargin)
  %
  % Given a fine polygon P and a simplified P_coarse, it shrinks P
  % with MCF until it no longer intersects P_coarse
  %
  % Input:
  %   P  M x 2 list of polygon vertex positions of the initial fine mesh
  %   P_coarse  m x 2 list of polygon vertex positions of the coarse mesh
  %   Optional:
  %     'delta_t' followed by timestep of the flow
  %     'flow_type' can be either 'curve' or 'surface' 
  %     or 'medial_attraction'
  % Output:
  %   P_all   Mx2xsteps lists of polygon vertex positions of the initial 
  %   fine mesh at each time step
  
  lambda = 1.0;
  flow_type = 'curve';
  % Parsing arguments
  ii = 1;
  while ii < numel(varargin)
      switch varargin{ii}
          case 'delta_t'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              lambda = varargin{ii};
          case 'flow_type'
              assert(ii+1<=numel(varargin));
              ii = ii+1;
              flow_type = varargin{ii};
          otherwise
              error('Unsupported parameter: %s',varargin{ii});
      end
      ii = ii+1;
  end
    
  % edges of the coarse mesh
  E = [(1:size(P,2))' [(2:size(P,2))';1]];
  
  % search for intersections between the coarse polygon edges and 
  % the fine polygon
  intersects = 0;
  for k=1:size(P_coarse,2)-1
      s = [P_coarse(:,k) P_coarse(:,k+1)];
      I = seg2poly_original(s,P);
      if ~isempty(I)
          intersects = 1;
          break;
      end
  end
  % check also for last segment
  s = [P_coarse(:,end) P_coarse(:,1)];
  I = seg2poly_original(s,P);
  if ~isempty(I)
      intersects = 1;
  end
  
  
  %initialize flow
  ii = 1;
  Sall(:,:,1) = P;
  
  % Create matrices and triangulayte depending on the flow
  switch flow_type
      case 'curve'
          A = adjacency_matrix(E);
          L = A-diag(sum(A,2));
      case 'surface'
          % first assert if there aren't duplicate vertices around each 
          % vertex
          for k=setdiff(-5:5,-1)
              min_dist_points = min(normrow([P(1,1:size(P,2))-P(1,mod(1+k:size(P,2)+k,size(P,2))+1);...
                  P(2,1:size(P,2))-P(2,mod(1+k:size(P,2)+k,size(P,2))+1)]'));
              assert(min_dist_points>1e-16,...
                  'ERROR: duplicate vertices found. Please provide other input');
          end
          [SV,SF] = triangle(P',E,[]);
          SV = [SV zeros(size(SV,1),1)];
          L = cotmatrix(SV,SF);
      case 'medial_attraction'
          
         A = adjacency_edge_cost_matrix(P',E);
         L = A-diag(sum(A,2));
         
         % number of explicit smoothnig steps
         k_max = 50;
         
         % implicit smoothing parameters
         s_L = 1.0;
         s_H = 0.005;
          
          % calculate medial axis of the coarse mesh
          pdeGeom = geomDataFromPolygon(P_coarse');
          [medialCurves,~,~] = computeMedialAxis(pdeGeom);
          hold on
          plotMedialCurves(medialCurves);
          nCurves = length(medialCurves);
          medialPoints = [];
          Emedial = [];
          controlPts = [];
          
          % sample medial axis and also genereate egdes
          for i = 1:nCurves
              curve = medialCurves(i);
              t = 0:(1/20):1;
              medialPoints = [medialPoints sampleBezier(curve,t)'];
              Emedial = [Emedial; ((i-1)*21+1:(i-1)*21+20)' ((i-1)*21+2:(i-1)*21+21)'];
              controlPts = [controlPts (curve.controlPoints)'];
          end
                    
          % for each point in the fine mesh, find its orthogonal
          % projection in the medial axis
%           P_medial_ortho = zeros(size(P,1),size(P,2));
%           for k=1:size(P,2)
%               
%               % P belongs to the medial axis, the projection is itself
% %               min(normrow((P(:,k)*ones(1,size(controlPts,2))-controlPts)'))
%               if min(normrow((P(:,k)*ones(1,size(controlPts,2))-controlPts)'))<1e-6
%                   P_medial_ortho(:,k) = P(:,k);
%                   continue;
%               end
%               
%               % Calculate normal
%               if k==1
%                   N = (([-(P(2,2)-P(2,1)) (P(1,2)-P(1,1))]/norm([-(P(2,2)-P(2,1)) (P(1,2)-P(1,1))])...
%                       +[-(P(2,1)-P(2,end)) (P(1,1)-P(1,end))]/norm([-(P(2,1)-P(2,end)) (P(1,1)-P(1,end))])))';
%                   N = N/norm(N);
%               elseif (k==size(P,2))
%                   N = (([-(P(2,end-1)-P(2,end)) (P(1,end-1)-P(1,end))]/norm([-(P(2,end-1)-P(2,end)) (P(1,end-1)-P(1,end))])...
%                       +[-(P(2,end)-P(2,1)) (P(1,end)-P(1,1))]/norm([-(P(2,end)-P(2,1)) (P(1,end)-P(1,1))])))';
%                   N = N/norm(N);
%               else
%                   N = (([-(P(2,k+1)-P(2,k)) (P(1,k+1)-P(1,k))]/norm([-(P(2,k+1)-P(2,k)) (P(1,k+1)-P(1,k))])...
%                       +[-(P(2,k)-P(2,k-1)) (P(1,k)-P(1,k-1))]/norm([-(P(2,k)-P(2,k-1)) (P(1,k)-P(1,k-1))])))';
%                   N = N/norm(N);
%               end
%               V0 = [medialPoints'; P(:,k)'];
%               V1 = [medialPoints'; P(:,k)'+100*N'];
%               [~,~,~,~,~,~,ColT] = collision_detection_2d_opt(V0,V1,[Emedial;size(V0,1) size(V0,1)]','last_only',1);
% 
%               if sum(sum(ColT))==0
%                   P_medial_ortho(:,k) = P(:,k);
%               else
%                   t_coeff = min(min(ColT+100000*(ColT==0)));
%                   P_medial_ortho(:,k) = (1-t_coeff)*P(:,k)+t_coeff*(P(:,k)+100*N);
%               end
%           end 
                  

  end
  
%   plot([P(1,:) P(1,1)],[P(2,:) P(2,1)],'o--k','LineWidth',3);
%   pcoarse = plot([P_coarse(1,:) P_coarse(1,1)],[P_coarse(2,:) P_coarse(2,1)],'o-r','LineWidth',3);
%   pc = plot([P(1,:) P(1,1)],[P(2,:) P(2,1)],'o-b','LineWidth',2);
%   input('')
%   delete(pc);
  
  while(intersects)
      
      ii = ii+1;
      
      switch flow_type
          case 'curve'
              A = adjacency_edge_cost_matrix(Sall(:,:,ii-1)',E);
              M = diag(sum(A,2)/2);
              M = M./max(diag(M));
              % cMCF
              Sall(:,:,ii) = ((M-lambda*L)\(M*Sall(:,:,ii-1)'))';
              P = Sall(:,:,ii);
          case 'surface'
              D = massmatrix(SV,SF,'full');
              S = D - lambda*L;
              SV = S\(D*SV);
              Sall(:,:,ii) = SV(1:size(P,2),1:2)';
              P = Sall(:,:,ii);
          case 'medial_attraction'
              % replace by the right flow
              
              % for each point in the fine mesh, find the closest in the medial
              % axis
              P_medial_closest = zeros(size(P,1),size(P,2));
              for k=1:size(P,2)
                  dist = 10000000;
                  for j=1:size(medialPoints,2)
                      if norm(P(:,k)-medialPoints(:,j))<dist
                          dist = norm(P(:,k)-medialPoints(:,j));
                          P_medial_closest(:,k) = medialPoints(:,j);
                      end
                  end
              end
              
%               A = adjacency_edge_cost_matrix(P',E);
%               L = A-diag(sum(A,2));
              
              P_inter_prev = P;
              
              t = min((ii-1)*lambda,1.0);
%               t = 1.0
              P_inter = (1-t)*P'+t*P_medial_closest';
%               % explicit smmothing
%               for k=1:k_max
%                   P_inter = P_inter+(t)*L*P_inter;
%               end
              %implicit smoothing
              P_inter = [s_L*L; s_H*speye(size(P_inter,1))]\[zeros(size(P_inter));s_H*P_inter];
              s_L = s_L;
              s_H = 1.05*s_H;
              
              Sall(:,:,ii) = P_inter';
              P = P_inter';
              
%               norm(P-P_inter_prev)
%               if norm(P-P_inter_prev)<1e-2
%                   disp('multiplying k_max by 2')
%                   k_max = round(2*k_max);
%               end
              
      end
      
      pc = plot([P(1,:) P(1,1)],[P(2,:) P(2,1)],'o-b','LineWidth',2);
      
      % check for intersections of the polygons
      intersects = 0;
      for k=1:size(P_coarse,2)-1
          s = [P_coarse(:,k) P_coarse(:,k+1)];
          I = seg2poly_original(s,P);
          % if points in I belong to P_coarse, ignore them
          j=1;
          while j<=size(I,2)
              if (norm(I(:,j)-P_coarse(:,k))<1e-6 || norm(I(:,j)-P_coarse(:,k+1))<1e-6)
                  I = [I(:,1:j-1) I(:,j+1:end)];
              else
                  j=j+1;
              end
          end
          if ~isempty(I)
              intersects = 1;
              break;
          end
      end
      
      % check also for last segment
      s = [P_coarse(:,end) P_coarse(:,1)];
      I = seg2poly_original(s,P);
      % if points in I belong to P_coarse, ignore them
      j=1;
      while j<=size(I,2)
          if (norm(I(:,j)-P_coarse(:,end))<1e-6 || norm(I(:,j)-P_coarse(:,1))<1e-6)
              I = [I(:,1:j-1) I(:,j+1:end)];
          else
              j=j+1;
          end
      end
      if ~isempty(I)
          intersects = 1;
      end
      
       
      ii
      input('');
%       pause(0.001);
      delete(pc);
      
  end
  
%   delete(pcoarse);
