function P_coarse_final = positions_expand_fine_pushing_coarse(P_all,P_coarse,delta_t,varargin)
% POSITIONS_EXPAND_FINE_PUSHING_COARSE
% positions_expand_fine_pushing_coarse(P_all,P_coarse,delta_t,varargin)
%
% Given the evolution of a fine polygon P (the last is assumed to be
% inside the coarse mesh), grow the fine polygon to its initial position
% pushing P_coarse in a way that results in no collisions.
%
% Input:
%   P_all  Mx2xsteps list of polygon vertex positions of the if the flown
%   fine polygon
%   P_coarse  m x 2 list of polygon vertex positions of the coarse mesh
%   delta_t   time step of the flow (will be the same time step for the
%   simulation
%   Optional:
%     'loop_collisions' followed by number of attempts to get to a
%     collision-free state at each time step of the physical simulation
%     'QPMethod' followed by either 'min_quad_with_fixed_active_set' or
%       'quadprog'
%     'energy' followed by either 'displacement_step',
%       'displacement_initial', 'dirichlet', 'laplacian2', 'surface_arap'
%     or 'proximity'
%     'simulation_steps' number of physical simulation steps to reach
%     one step back of the flow
%     'back_factor' factor that determines how much the vertices
%     of the coarse mesh will try to go back to their initial positions
%     when 'energy' is 'displacement_initial'
% Output:
%   P_coarse_final   m x 2 list of polygon vertex positions of the
%   resulting cage

loop_collisions = 100;
qp_method = 'quadprog';
energy = 'displacement_step';
simulation_steps = 1;
back_factor = 1.0;
% Parsing arguments
ii = 1;
while ii < numel(varargin)
    switch varargin{ii}
        case 'loop_collisions'
            assert(ii+1<=numel(varargin));
            ii = ii+1;
            loop_collisions = varargin{ii};
        case 'QPMethod'
            ii = ii+1;
            assert(ii<=numel(varargin));
            qp_method = varargin{ii};
        case 'energy'
            ii = ii+1;
            assert(ii<=numel(varargin));
            energy = varargin{ii};
        case 'simulation_steps'
            ii = ii+1;
            assert(ii<=numel(varargin));
            simulation_steps = varargin{ii};
        case 'back_factor'
            ii = ii+1;
            assert(ii<=numel(varargin));
            back_factor = varargin{ii};
        otherwise
            error('Unsupported parameter: %s',varargin{ii});
    end
    ii = ii+1;
end

% fine mesh
OV = P_all(:,:,1)';
E = [1:size(OV,1);2:size(OV,1) 1]';
% coarse  mesh
CV = P_coarse';
CE = [1:size(CV,1);2:size(CV,1) 1]';
% Shrunk fine mesh
V = P_all(:,:,end)';

pov = plot(reshape(OV(E',1),2,[]),reshape(OV(E',2),2,[]),'o--k','LineWidth',1);
hold on;
pv = plot(reshape(V(E',1),2,[]),reshape(V(E',2),2,[]),'o-r','LineWidth',3);
pc = plot(reshape(CV(CE',1),2,[]),reshape(CV(CE',2),2,[]),'o-b','LineWidth',3);
title(sprintf('Get ready! '));
hold off;
axis equal;
input('Ready. Set. ');
fprintf('Go!\n');

t_step = delta_t;
ts = fliplr(size(P_all,3)*delta_t:-delta_t:0);
ts = ts(ts>0);

% initial coarse mesh is used during the simulation
CV_init = CV;

% edges for all vertices of both polygons
E_all = [E;CE+size(V,1)];

pcol = [];
ncol = [];

% number of steps
k = size(P_all,3);
t_step_simulation = 1/simulation_steps;

% define mass and cotangent matrices w.r.t. initial embedding, for both
% fine (constrained) and coarse (unconstrained) meshes
Adjc = adjacency_edge_cost_matrix(V,E);
Mc_in = diag(sum(Adjc,2)/2);
Mc_in = Mc_in./max(diag(Mc_in));
Mc_in = blkdiag(Mc_in,Mc_in);
Adjun = adjacency_edge_cost_matrix(CV,CE);
Mun_in = diag(sum(Adjun,2)/2);
Mun_in = Mun_in./max(diag(Mun_in));
Mun_in = blkdiag(Mun_in,Mun_in);
Lc_in = Adjc-diag(sum(Adjc,2));
Lc_in = blkdiag(Lc_in,Lc_in);
Lun_in = Adjun-diag(sum(Adjun,2));
Lun_in = blkdiag(Lun_in,Lun_in);


for t = ts(1:end-1)
    
    CV_loop = [CV(:,1);CV(:,2)];
    
    for k_simulation=1:simulation_steps
        
        % previous configuration of all vertices (for collision detection)
        V_all_prev = [V;CV];
        
        % update positions of the fine (constrained) mesh
        V = (P_all(:,:,k-1)'-P_all(:,:,k)')*t_step_simulation + V;

        % define mass and cotangent matrices for both meshes
        Adjc = adjacency_edge_cost_matrix(V,E);
        Mc = diag(sum(Adjc,2)/2);
        Mc = Mc./max(diag(Mc));
        Mc = blkdiag(Mc,Mc);
        Adjun = adjacency_edge_cost_matrix(CV,CE);
        Mun = diag(sum(Adjun,2)/2);
        Mun = Mun./max(diag(Mun));
        Mun = blkdiag(Mun,Mun);
        Lc = Adjc-diag(sum(Adjc,2));
        Lc = blkdiag(Lc,Lc);
        Lun = Adjun-diag(sum(Adjun,2));
        Lun = blkdiag(Lun,Lun);
        
        % System Matrix
        switch energy
            case 'displacement_step'
                A = Mun;
            case 'displacement_initial'
                A = Mun_in;
            case 'dirichlet'
                % Dirichlet energy of the displacement
                A = Mun*Lun;
            case 'dirichlet_initial'
                % Dirichlet energy of the displacement w.r.t. initial
                % embedding
                A = Mun_in*Lun_in;
        end
        
        % rhs
        switch energy
            case 'displacement_step'
                b = -Mun*[CV(:,1);CV(:,2)];
            case 'displacement_initial'
                b = -Mun_in*[CV_init(:,1);CV_init(:,2)];
            case 'dirichlet'
                b = -Mun*Lun*[CV(:,1);CV(:,2)];
            case 'dirichlet_initial'
                b = -Mun_in*Lun_in*[CV_init(:,1);CV_init(:,2)];
        end
        
        % Calculate new positions
        qpA = A;
        qpB = b;
        qpknown = [];
        qpY = [];
        qpAeq = [];
        qpBeq = [];
        qpAieq = [];
        qpBieq = [];
        
        % To plot progress:
        CV_loop_prev = CV_loop;
        switch qp_method
            case 'min_quad_with_fixed_active_set'
                CV_loop = min_quad_with_fixed_active_set(qpA,qpB,qpknown,qpY,qpAeq,qpBeq,qpAieq,qpBieq);
            case 'quadprog'
                CV_loop = quadprog( ...
                    qpA,qpB,qpAieq,qpBieq, ...
                    [qpAeq;sparse(1:numel(qpknown),qpknown,1,numel(qpknown),size(qpA,2))], ...
                    [qpBeq;qpY],[],[]);
        end

        
         % update list of all vertex positions (for collision detection)
        V_all = [V;[CV_loop(1:end/2,:) CV_loop(end/2+1:end,:) ]];
        
        % initialize contraint matrix and rhs as empty
        NB = [];
        rhs = [];

        % find collisions
        [Collisions,ColX,ColY,ColNX,ColNY,ColU] = collision_detection_2d_opt(V_all_prev,V_all,E_all');
        % ignoring intersections in the fine mesh
        Collisions(1:size(V,1),1:size(E,1)) = sparse(size(V,1),size(E,1));
        
        % constraint matrices Nun and Nc
        [N,B,Acon] = positions_response_constraint_matrix(E_all,Collisions,ColNX,ColNY,ColU,V,E);
        % add constraints
        NB = [NB;N*B];
        rhs = [rhs; N*Acon*[V(:,1);V(:,2)]];
        
        % x and y coordinates of the colliding points (only for plotting)
        colliding_points = full([ColX(find(Collisions)) ColY(find(Collisions))]);
        % collising normals (just for plotting
        colliding_normals = full([ColNX(find(Collisions)) ColNY(find(Collisions))]);
            
        loop_iter = 0;        
        % plot progress of the coarse mesh
        progress = norm(CV_loop_prev-CV_loop,inf);
        
        hold on;
        delete(pv);
        pv = plot(reshape(V(E',1),2,[]),reshape(V(E',2),2,[]),'o-r','LineWidth',3);
        delete(pc);
        CV_loop_drawing = [CV_loop(1:end/2) CV_loop(end/2+1:end)];
        pc = plot(reshape(CV_loop_drawing(CE',1),2,[]),reshape(CV_loop_drawing(CE',2),2,[]),'o-b','LineWidth',3);
        % plot collisions
        delete(pcol);
        pcol = plot(colliding_points(1:end/2),colliding_points(end/2+1:end),'*g');
        % plot normals
        delete(ncol);
        ncol = vectarrow_new([colliding_points(1:end/2)' colliding_points(end/2+1:end)'],...
            0.1*[colliding_normals(1:end/2)' colliding_normals(end/2+1:end)']);
        set(ncol,'Color','g','Linewidth',1.5);
        title(sprintf('collision attempt: %d/%d, #constraints = %d, vel. change = %.6f \n simulation step: %d/%d, flow step: %d/%d, energy:%s ', ...
            loop_iter, loop_collisions, size(rhs,1), progress, k_simulation,simulation_steps,...
            size(P_all,3)-(k-1), size(P_all,3)-1, strrep(energy, '_', '\_')));
        drawnow;
%         pause(0.005);
%         input('');
        hold off;

        % - repeat until feasible or the coarse mesh stops progressing
        progress = 1.0;
        while (sum(sum(Collisions))>0 && loop_iter<loop_collisions && progress >= 0.0)

            
            qpA = A;
            qpB = b;
            qpknown = [];
            qpY = [];
            qpAeq = [];
            qpBeq = [];
            % Matlab inequalities are less or equal (so flip signs)
            qpAieq = -NB;
            qpBieq = -rhs-1e-2;

            
            % To plot progress:
            CV_loop_prev = CV_loop;
            switch qp_method
                case 'min_quad_with_fixed_active_set'
                    CV_loop = min_quad_with_fixed_active_set(qpA,qpB,qpknown,qpY,qpAeq,qpBeq,qpAieq,qpBieq);
                case 'quadprog'
                    CV_loop = quadprog( ...
                        qpA,qpB,qpAieq,qpBieq, ...
                        [qpAeq;sparse(1:numel(qpknown),qpknown,1,numel(qpknown),size(qpA,2))], ...
                        [qpBeq;qpY],[],[]);
            end
            
            % update list of all vertex positions (for collision detection)
            V_all = [V;[CV_loop(1:end/2,:) CV_loop(end/2+1:end,:)]];
                
%             disp('constraint satisfaction')
%             qpAieq*Deltavel-qpBieq
            
            % find collisions
            [Collisions,ColX,ColY,ColNX,ColNY,ColU] = collision_detection_2d_opt(V_all_prev,V_all,E_all');
            % ignoring intersections in the fine mesh
            Collisions(1:size(V,1),1:size(E,1)) = sparse(size(V,1),size(E,1));

            % constraint matrices Nun and Nc
            [N,B,Acon] = positions_response_constraint_matrix(E_all,Collisions,ColNX,ColNY,ColU,V,E);
            % add constraints
            NB = [NB;N*B];
            rhs = [rhs; N*Acon*[V(:,1);V(:,2)]];
            loop_iter = loop_iter+1;

            % x and y coordinates of the colliding points (only for plotting)
            colliding_points = full([ColX(find(Collisions)) ColY(find(Collisions))]);
            % collising normals (just for plotting
            colliding_normals = full([ColNX(find(Collisions)) ColNY(find(Collisions))]);

            
            % plot progress of the coarse mesh
            progress = norm(CV_loop_prev-CV_loop,inf);

            hold on;
            delete(pv);
            pv = plot(reshape(V(E',1),2,[]),reshape(V(E',2),2,[]),'o-r','LineWidth',3);
            delete(pc);
            CV_loop_drawing = [CV_loop(1:end/2) CV_loop(end/2+1:end)];
            pc = plot(reshape(CV_loop_drawing(CE',1),2,[]),reshape(CV_loop_drawing(CE',2),2,[]),'o-b','LineWidth',3);
            % plot collisions
            delete(pcol);
            pcol = plot(colliding_points(1:end/2),colliding_points(end/2+1:end),'*g');
            % plot normals
            delete(ncol);
            ncol = vectarrow_new([colliding_points(1:end/2)' colliding_points(end/2+1:end)'],...
                0.1*[colliding_normals(1:end/2)' colliding_normals(end/2+1:end)']);
            set(ncol,'Color','g','Linewidth',1.5);
            title(sprintf('collision attempt: %d/%d, #constraints = %d, vel. change = %.6f \n simulation step: %d/%d, flow step: %d/%d, energy: %4s', ...
                loop_iter, loop_collisions, size(rhs,1), progress, k_simulation,simulation_steps,...
                size(P_all,3)-(k-1), size(P_all,3)-1,strrep(energy, '_', '\_')));
            drawnow;
%             pause(0.005);
%             input('');
            hold off;

        end
        
        CV = [CV_loop(1:end/2,:) CV_loop(end/2+1:end,:)];
         
    end
    
    k = k - 1;
    
end

P_coarse_final = CV';