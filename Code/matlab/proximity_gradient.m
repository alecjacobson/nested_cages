function grad = proximity_gradient(P,P_coarse,u,inds)
% PROXIMITY_GRADIENT
% H = proximity_gradient(P,P_coarse,u,inds)
%
% Calculate the hessian to minimize proximity energy
%
% Input:
%   P  2xn_all list of polygon vertex positions of the fine mesh
%   P_coarse 2xn_coarse list of polygon vertex positions of the fine mesh
%   u:  1xn_all list orthogonal projection coefficients for each point in P
%   inds: indices of the points in P_coarse in P's order
% Output:
%   grad: (n_all+n+coarse)x1 gradient vector

P = P';
P_coarse = P_coarse';
grad = zeros(2*(size(P_coarse,1)+size(P,1)),1);

for k=1:size(inds,2)-1
    
    for j=inds(k):inds(k+1)-1
        
        grad(k) = 2*(1-u(j))^2*P_coarse(k,1) + 2*(1-u(j))*u(j)*P_coarse(k+1,1)-2*u(j)*P(j,1);
        grad(k+size(P_coarse,1)+size(P,1)) = 2*(1-u(j))^2*P_coarse(k,2) + 2*(1-u(j))*u(j)*P_coarse(k+1,2)-2*u(j)*P(j,2);
        
        grad(k+1) = 2*u(j)^2*P_coarse(k+1,1) + 2*(1-u(j))*u(j)*P_coarse(k,1)-2*(1-u(j))*P(j,1);
        grad(k+1+size(P_coarse,1)+size(P,1)) = 2*u(j)^2*P_coarse(k+1,2) + 2*(1-u(j))*u(j)*P_coarse(k,2)-2*(1-u(j))*P(j,2);
        
    end
    
end

grad = zeros-grad;