function V_coarse_final = final_velocityfilter_step_project(V0,F0,V_coarse,F_coarse,varargin)
% To-do: add the other energies (currently we only have volume)
% FINAL_VELOCITYFILTER_STEP_PROJECT
% V_coarse_final = final_velocityfilter_step_project(V0,F0,V_coarse,F_coarse,varargin)
%
% Given (V0,F0) and an initial coarse cage (V_coarse,F_coarse) (i.e., the
% coarse mesh is already fesible as a cage), step (V_coarse,F_coarse)
% towards the gradient of some energy and project it back to the feasible 
% set with velocityfilter
%
% Input:
%   V0  (#vertices)x3 list of mesh vertex positions of the initial fine mesh
%   F0  (#faces)x3 list of vertex indices that form each face of the
%   initial mesh
%   V_coarse   (#vertices_cage)x3 list of mesh vertex positions of the 
%   coarse mesh
%   F_coarse   (#faces_cage)x3 list of vertex indices that form each face
%   of the coarse mesh
%   Optional:
%     'energy': followed by either 'displacement_step' or
%       'displacement_initial' or 'symmetry_x' or 
%       'volume' (default)
%     'beta': step size towards gradient direction
% Output:
%   V_coarse_final (#vertices_cage)x3 list of mesh vertex positions of the
%   resulting cage

IF = intersect_other(V0,F0,V_coarse,F_coarse,'FirstOnly',true);
if size(IF,1)>0
    error('meshes intersect. Aborting');
end

energy = 'volume';
beta = 1e-3;
% parsing arguments
ii = 1;
while ii < numel(varargin)
    switch varargin{ii}
        case 'energy'
            ii = ii+1;
            assert(ii<=numel(varargin));
            energy = varargin{ii};
        case 'beta'
            ii = ii+1;
            assert(ii<=numel(varargin));
            beta = varargin{ii};
        otherwise
            error('Unsupported parameter: %s',varargin{ii});
    end
    ii = ii+1;
end
disp('started finalvelocityfilter_step_project');
V = V0;
F_all = [F0; size(V0,1)+F_coarse];
beta_init = beta;

if (strcmp(energy,'volume'))
                
            % configuration is already feasible
            CV = V_coarse;
            CF = F_coarse;
            [~,min_energy] = centroid(CV,CF);
            CV_opt = CV;
            fprintf('initial energy = %g\n', min_energy);

            while(beta>1e-6)
                
                V_all_prev = [V;CV];
                
                areas = doublearea(CV,CF)/2;
                N = normalizerow(normals(CV,CF));

                grad_vol = zeros(size(CV));
                for i=1:size(CV,1)
                    face_idx = mod(find(CF==i)-1,size(CF,1))+1;
                    grad_vol(i,:) = sum([areas(face_idx) areas(face_idx) areas(face_idx)].*N(face_idx,:));
                end
                
                fprintf('number of known vertices: %d\n',size(V,1));
                % step+project
                Z = velocity_filter_mex(V_all_prev,[V;CV-beta*grad_vol],F_all,size(V,1));
                gamma = 1;
                CV_attempt = CV+gamma*(Z(size(V,1)+1:end,:)-CV);
                
                [~,cur_energy] = centroid(CV_attempt,CF);
                
                % if energy decreased, 'save' this state
                if (cur_energy<min_energy)
                    % update meshes
                    CV = CV_attempt;
                    % save optimal result to output
                    CV_opt = CV;
                    min_energy = cur_energy;
                    beta = beta_init;
                    fprintf('energy decreased: energy = %g\n', cur_energy);
                else
                    beta = 0.5*beta;
                    fprintf('energy increased to energy = %g. Changing beta = %g\n', cur_energy, beta);
                end
                
            end
            
            CV = CV_opt;
                
end

V_coarse_final = CV;