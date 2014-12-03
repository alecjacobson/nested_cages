% fine mesh
nf = 20;
theta = linspace(0,2*pi,nf+1)';
theta = theta(1:end-1);
OV = [cos(theta) sin(theta)];
E = [1:size(OV,1);2:size(OV,1) 1]';
% coarse  mesh
nc = 6;
theta = linspace(0,2*pi,nc+1)';
theta = theta(1:end-1);
CV = [cos(theta) sin(theta)];
CE = [1:size(CV,1);2:size(CV,1) 1]';

% Shrink fine mesh
r = 0.9*min(sqrt(sum(barycenter(CV,CE).^2,2)));
V = r*OV;
      
pov = plot(reshape(OV(E',1),2,[]),reshape(OV(E',2),2,[]),'o--k','LineWidth',1);
hold on;
  pv = plot(reshape(V(E',1),2,[]),reshape(V(E',2),2,[]),'o-r','LineWidth',3);
  pc = plot(reshape(CV(CE',1),2,[]),reshape(CV(CE',2),2,[]),'o-b','LineWidth',3);
hold off;
axis equal;

t_step = 0.1;
ts = fliplr(1:-t_step:0);
ts = ts(ts>0);
% constant velocity
Vvel = OV-V;
CVvel = zeros(size(CV,1),2);

% edges for all vertices of both polygons
E_all = [E;CE+size(V,1)];
pcol = [];

for t = ts
  fprintf('t: %g\n',t);
  % previous configuration of all vertices
  V_all_prev = [V;CV];
  % Fine mesh progresses regardless: Vvel is fixed
  V = Vvel*t_step + V;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TODO: Solve for new CV and new CVvel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % - solve for unconstrained update of CVvel

  % Set system matrix and rhs: eq. (3) from Otaduy et al. 09.
  % Mass Matrix
  Adj = adjacency_edge_cost_matrix([V;CV],[E; CE+size(V,1)]);
  M = diag(sum(Adj,2)/2);
  M = M./max(diag(M));
  M = blkdiag(M,M);
  % System Matrix (F = -2*(x-x0), so derivative is -2*eye)
  A = M - (t_step^2)*(-2*speye(2*size([V;CV],1)));
  % rhs: F(x0) = 0, so the first term is canceled
  b = M*[Vvel(:,1);CVvel(:,1);Vvel(:,2);CVvel(:,2)];
  % initialize matrix J and vector g0 in (5) from Otaduy et al. 09 as empty
  J = [];
  g0 = [];
  
  % Calculate velocity update 
  Deltavel = min_quad_with_fixed_active_set(A,-b,...
      [1:size(Vvel,1) (size(Vvel,1)+size(CVvel,1))+1:(size(Vvel,1)+size(CVvel,1))+size(Vvel,1)],...
      [Vvel(1:size(Vvel,1),1);Vvel(1:size(Vvel,1),2)]);
  % Select velocities from vertices that should move
  DeltaCVvel = [Deltavel(size(V,1)+1:size(V,1)+size(CV,1))...
      Deltavel((size(V,1)+size(CV,1))+size(V,1)+1:(size(V,1)+size(CV,1))+size(V,1)+size(CV,1))];
  % update velocitiy
  CVvel = CVvel + DeltaCVvel;
  % Calculate position update
  DeltaCV = t_step*(CVvel);
  % update positions
  CV = CV + DeltaCV;
  
  % - TODO: tolerance-based collision detection from CV to initialize constraint set
  
  % update list of all vertex velocities
  V_all = [V;CV];
  % find collisions
  [Collisions,ColX,ColY,ColNX,ColNY,ColU] = collision_detection_2d_opt(V_all_prev,V_all,E_all');
  % constraint matrix J
  [J_new,Dg_Dp,Dp_Dv] = response_constraint_matrix(E_all,Collisions,ColNX,ColNY,ColU);
  % rhs vector from inequality (5). It includes -1/delta_t
  rhs = response_constraint_rhs(E_all,Collisions,ColNX,ColNY,ColU,V_all_prev,t_step);
  % add constraints
  J = [J;J_new];
  g0 = [g0; rhs];
  
%   % test equation (5)
%   % all velocities
%   Vvel_all = [Vvel; CVvel]
%   Vvel_all_vec = [Vvel_all(:,1);Vvel_all(:,2)];
%   % obs.: velocities at the beginning of the timestep are zero
%   full([J*zeros(size(Vvel_all_vec)) rhs J*Vvel_all_vec rhs])
  
  % x and y coordinates of the colliding points (only for plotting)
  colliding_points = [ColX(find(Collisions)) ColY(find(Collisions))];
  
  % - repeat until feasible
  while (sum(sum(Collisions))>0)
      disp('there are collisions');
  %   - Collision response via MLCP: change in CVvel
        Vvel_all = [Vvel(:,1); CVvel(:,1); Vvel(:,2); CVvel(:,2)];
        Deltavel = min_quad_with_fixed_active_set(A,sparse(size(Vvel_all,1),1),...
            [1:size(Vvel,1) (size(Vvel,1)+size(CVvel,1))+1:(size(Vvel,1)+size(CVvel,1))+size(Vvel,1)],...
            [Vvel(1:size(Vvel,1),1);Vvel(1:size(Vvel,1),2)],[],[],-J,-(g0-J*Vvel_all));
        % Select velocities from vertices that should move
        DeltaCVvel = [Deltavel(size(V,1)+1:size(V,1)+size(CV,1))...
            Deltavel((size(V,1)+size(CV,1))+size(V,1)+1:(size(V,1)+size(CV,1))+size(V,1)+size(CV,1))];
  %   - Compute tentative update for CVvel and CV
        % update velocitiy
        CVvel = CVvel + DeltaCVvel;
        % Calculate position update
        DeltaCV = t_step*(CVvel);
        % update positions
        CV = CV + DeltaCV;
  %   - Collision detection update constraints
        % update list of all vertex velocities
        V_all = [V;CV];
        % find collisions
        [Collisions,ColX,ColY,ColNX,ColNY,ColU] = collision_detection_2d_opt(V_all_prev,V_all,E_all');
        % constraint matrix J
        [J_new,Dg_Dp,Dp_Dv] = response_constraint_matrix(E_all,Collisions,ColNX,ColNY,ColU);
        % rhs vector from inequality (5). It includes -1/delta_t
        rhs = response_constraint_rhs(E_all,Collisions,ColNX,ColNY,ColU,V_all_prev,t_step);
        % add constraints
        J = [J;J_new];
        g0 = [g0; rhs];
  end
  
  % reset constraints to none and CVvel to zero
  J = [];
  g0 = [];
  CVvel = zeros(size(CV,1),2);
  
  hold on;
    delete(pv);
    pv = plot(reshape(V(E',1),2,[]),reshape(V(E',2),2,[]),'o-r','LineWidth',3);
    delete(pc);
    pc = plot(reshape(CV(CE',1),2,[]),reshape(CV(CE',2),2,[]),'o-b','LineWidth',3);
    % plot collisions
    delete(pcol);
    pcol = plot(colliding_points(1:end/2),colliding_points(end/2+1:end),'*g');
    drawnow;
    %pause(0.5);
  hold off;
  input('');
end

