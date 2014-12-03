function [u,inds] = ortho_proj_coefficients(P,P_coarse)
% ORTHO_PROJ_COEFFICIENTS
% [u,inds] = ortho_proj_coefficients(P,P_coarse)
%
% Identify points in P_coarse that belong to P and for each 
% point in P find the orthogonal projection coefficent
% in the corresponding segment of P_coarse
%
% Input:
%   P  2xn_all list of polygon vertex positions of the fine mesh
%   P_coarse 2xn_coarse list of polygon vertex positions of the fine mesh
% Output:
%   u:  1xn_all list orthogonal projection coefficients for each point in P
%   inds: indices of the points in P_coarse in P's order

% transpose
P = P';
P_coarse = P_coarse';

% number of vertices for both meshes
n_all = size(P,1);
n_coarse = size(P_coarse,1);

% find points in P that belong to P_coarse
inds = [];
for i=1:n_coarse
    for j=1:n_all
        if (norm(P_coarse(i,:)-P(j,:))<1e-16)
            inds = [inds j];
            break;
        end
    end
end

% for each point in P, find the coefficient of its orthogonal projection
% on the corresponding segment of P
u = [];
for i=1:n_coarse-1
    k1 = inds(i);
    k2 = inds(i+1);
    for k=k1:k2-1
        v = P(k,:)-P_coarse(i,:);
        v_coarse = P_coarse(i+1,:)-P_coarse(i,:);
        u = [u dot(v,v_coarse)/(norm(v_coarse)^2)];
    end
end
