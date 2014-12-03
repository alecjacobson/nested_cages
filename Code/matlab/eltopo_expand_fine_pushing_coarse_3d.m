function [V_coarse_final,V_new] = eltopo_expand_fine_pushing_coarse_3d(P_all,F0,V_coarse,F_coarse,delta_t,varargin)
% ELTOPO_EXPAND_FINE_PUSHING_COARSE_3D
% eltopo_expand_fine_pushing_coarse_3d(P_all,F0,V_coarse,F_coarse,delta_t,varargin)
%
% Given the evolution of a fine mesh V (the last is assumed to be
% inside the coarse mesh), grow the fine mesh to its initial position
% pushing (V_coarse,F_coarse) in a way that results in no collisions,
% using the ElTopo library.
%
% Input:
%   P_all  (#vertices)x3xsteps list of mesh vertex positions of the initial fine mesh
%   F0  (#faces)x3 list of vertex indices that form each face of the
%   initial mesh
%   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
%   coarse mesh
%   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
%   of the coarse mesh
%   delta_t   time step of the flow (will be the same time step for the
%   simulation
%   Optional:
%     'simulation_steps' number of physical simulation steps to reach
%     one step back of the flow
%     'last_simulation' multiple to divide the last time step
%     where all levels are pushed simultaneously
% Output:
%   V_coarse_final (#vertices_cage)x3 list of mesh vertex positions of the
%   resulting cage

simulation_steps = 1;
last_simulation = 0;
energy = 'displacement_step';
% parsing arguments
ii = 1;
while ii < numel(varargin)
    switch varargin{ii}
        case 'simulation_steps'
            ii = ii+1;
            assert(ii<=numel(varargin));
            simulation_steps = varargin{ii};
        case 'last_simulation'
            ii = ii+1;
            assert(ii<=numel(varargin));
            last_simulation = varargin{ii};
        case 'energy'
            ii = ii+1;
            assert(ii<=numel(varargin));
            energy = varargin{ii};
        otherwise
            error('Unsupported parameter: %s',varargin{ii});
    end
    ii = ii+1;
end

% fine mesh
F = F0;
% coarse  mesh
CV = V_coarse;
CF = F_coarse;
% Shrunk fine mesh
V = P_all(:,:,end);

axis equal;
pc = [];
pv = [];
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

input('Ready. Set. ');
fprintf('Go!\n');

cla;

ts = fliplr(size(P_all,3)*delta_t:-delta_t:0);
ts = ts(ts>0);
% number of steps
k = size(P_all,3);
t_step_simulation = 1/simulation_steps;
% faces for all vertices of both meshes (for collision detection)
F_all = [F;CF+size(V,1)];

% define initial cage
CV_0 = CV;

eps_proximity = 1e-4;

for t = ts(1:end-1)
    
    if (k==2 && last_simulation)
        
        simulation_steps = simulation_steps*last_simulation;
        t_step_simulation = 1/simulation_steps;

%          eps_proximity = eps_proximity*last_simulation;
        
    end
    
    for k_simulation=1:simulation_steps
                
        % previous configuration of all vertices (for collision detection)
        V_all_prev = [V;CV];    
        
        % update positions of the fine (constrained) mesh
%         V = (P_all(:,:,k-1)-P_all(:,:,k))*t_step_simulation + V;
        V = (P_all(:,:,k-1)-P_all(:,:,k))*(k_simulation*t_step_simulation) + P_all(:,:,k); 
        
        if strcmp(energy,'displacement_step')
            
            % run ElTopo
            writeOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/meshes/V0_sphere_162.obj',V_all_prev,F_all);
            writeOBJ('/Users/Leo/PHD_Work/Cage_Generation_2013/code/eltopo/eltopo3d/meshes/V1_sphere_162.obj',[V;CV],F_all);
            V_eltopo = collide_eltopo_mex(V_all_prev,F_all,[V;CV],size(V,1),eps_proximity);
            CV = V_eltopo(size(V,1)+1:end,:);
        
            disp('max difference between ElTopo and constrained mesh')
            max(abs(V_eltopo(1:size(V,1),:)-V))
            
        elseif strcmp(energy,'displacement_initial')
            
            loop_iter = 1;
            V_eltopo = V_all_prev;
            tol = 5e-3;
            back_factor = 1.0;
            
            while((loop_iter==1 || max(normrow(V_all_prev-V_eltopo))>tol) && (back_factor>=0.0))
                
                V_all_prev = V_eltopo;
                CV_prev = V_eltopo(size(V,1)+1:end,:);
%                 CV_interp = back_factor*(CV_0-CV_prev)*(k_simulation*t_step_simulation) + CV_prev; 
                CV_interp = CV_0;
                
                % plot postions attempted
                if (loop_iter==1)
                    hold on;
                    delete(pc);
                    pc = trisurf(CF,CV_interp(1:end,1),CV_interp(1:end,2),CV_interp(1:end,3),'FaceColor',[0.0 0.5 0.0],'FaceAlpha',0.1);
                    delete(pv);
                    pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
                    title(sprintf('simulation step: %d/%d, flow step: %d/%d', ...
                        k_simulation,simulation_steps,...
                        size(P_all,3)-(k-1), size(P_all,3)-1));
                    drawnow;
                    input('');
                end
                
                V_eltopo = collide_eltopo_mex(V_all_prev,F_all,[V;CV_interp],size(V,1),eps_proximity);
                
                loop_iter = loop_iter+1;
                
                V = V_eltopo(1:size(V,1),:);
                CV = V_eltopo(size(V,1)+1:end,:);
                
                % plot partial results
                hold on;
                delete(pc);
                pc = trisurf(CF,CV(1:end,1),CV(1:end,2),CV(1:end,3),'FaceColor',[0.0 0.0 0.0],'FaceAlpha',0.1);
                delete(pv);
                pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
                title(sprintf('simulation step: %d/%d, flow step: %d/%d', ...
                    k_simulation,simulation_steps,...
                    size(P_all,3)-(k-1), size(P_all,3)-1));
                
%                 if (back_factor==0.0)
%                     break;
%                 end
                
%                 if mod(loop_iter,10)==0
%                     back_factor = back_factor-0.05
%                 end
                
%                 disp('progress')
%                 max(normrow(V_all_prev-V_eltopo))
                
                drawnow;
                input('');
                
            end
            
        end
        
        V = V_eltopo(1:size(V,1),:);
        
        % plot result of unconstrained minimization
        hold on;
        delete(pc);
        pc = trisurf(CF,CV(1:end,1),CV(1:end,2),CV(1:end,3),'FaceColor',[0.5 0.0 0.0],'FaceAlpha',0.1);
        delete(pv);
        pv = trisurf(F,V(:,1,end),V(:,2,end),V(:,3,end),'FaceColor',[0.0 0.0 0.8],'FaceAlpha',0.2);
        title(sprintf('simulation step: %d/%d, flow step: %d/%d', ...
             k_simulation,simulation_steps,...
            size(P_all,3)-(k-1), size(P_all,3)-1));
        drawnow;
        input('');
        hold off;
        
    end
    
    k = k-1;
    
end

V_new = V;
V_coarse_final = CV;